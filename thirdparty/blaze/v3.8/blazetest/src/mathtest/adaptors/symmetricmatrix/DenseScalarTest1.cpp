//=================================================================================================
/*!
//  \file src/mathtest/adaptors/symmetricmatrix/DenseScalarTest1.cpp
//  \brief Source file for the SymmetricMatrix dense scalar test (part 1)
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
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/util/policies/ArrayDelete.h>
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
/*!\brief Test of the SymmetricMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the SymmetricMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseScalarTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor (StaticMatrix)
   {
      test_ = "Row-major SymmetricMatrix default constructor (StaticMatrix)";

      const blaze::SymmetricMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > sym;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 0UL );
   }

   // Default constructor (HybridMatrix)
   {
      test_ = "Row-major SymmetricMatrix default constructor (HybridMatrix)";

      const blaze::SymmetricMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::rowMajor> > sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Default constructor (DynamicMatrix)
   {
      test_ = "Row-major SymmetricMatrix default constructor (DynamicMatrix)";

      const ST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   // Size constructor (HybridMatrix)
   {
      test_ = "Row-major SymmetricMatrix size constructor (HybridMatrix)";

      const blaze::SymmetricMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::rowMajor> > sym( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 0UL );
   }

   // Size constructor (DynamicMatrix)
   {
      test_ = "Row-major SymmetricMatrix size constructor (DynamicMatrix)";

      const ST sym( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Row-major list initialization
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Row-major SymmetricMatrix initializer list constructor (complete list)";

      const ST sym{ { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Row-major SymmetricMatrix initializer list constructor (incomplete list)";

      const ST sym{ { 1, 2, 3 }, { 2, 4 }, { 3, 0, 6 } };

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array initialization
   //=====================================================================================

   // Dynamic array initialization constructor
   {
      test_ = "Row-major SymmetricMatrix dynamic array initialization constructor";

      std::unique_ptr<int[]> array( new int[9] );
      array[0] = 1;
      array[1] = 2;
      array[2] = 3;
      array[3] = 2;
      array[4] = 4;
      array[5] = 0;
      array[6] = 3;
      array[7] = 0;
      array[8] = 6;
      const ST sym( 3UL, array.get() );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Static array initialization constructor
   {
      test_ = "Row-major SymmetricMatrix static array initialization constructor";

      const int array[3][3] = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };
      const ST sym( array );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major custom matrix constructors
   //=====================================================================================

   // Custom matrix constructor (ElementType*, size_t)
   {
      test_ = "Row-major SymmetricMatrix custom matrix constructor (ElementType*, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[5UL] );
      memory[1] = 1;
      memory[2] = 2;
      memory[3] = 2;
      memory[4] = 1;
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( memory.get()+1UL, 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 ||
          sym(1,0) != 2 || sym(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, size_t)
   {
      test_ = "Row-major SymmetricMatrix custom matrix constructor (ElementType*, size_t, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[11UL] );
      memory[1] = 1;
      memory[2] = 2;
      memory[6] = 2;
      memory[7] = 1;
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( memory.get()+1UL, 2UL, 5UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 ||
          sym(1,0) != 2 || sym(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Row-major SymmetricMatrix copy constructor (0x0)";

      const ST sym1;
      const ST sym2( sym1 );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Row-major SymmetricMatrix copy constructor (3x3)";

      ST sym1( 3UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      const ST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Row-major SymmetricMatrix move constructor (0x0)";

      ST sym1;
      ST sym2( std::move( sym1 ) );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Row-major SymmetricMatrix move constructor (3x3)";

      ST sym1( 3UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      ST sym2( std::move( sym1 ) );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Row-major SymmetricMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;
      const ST sym( mat );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Conversion constructor (symmetric)
   {
      test_ = "Row-major SymmetricMatrix conversion constructor (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    {  7,  0, 3 } } );

      const ST sym( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 0 ||
          sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-symmetric)
   {
      test_ = "Row-major SymmetricMatrix conversion constructor (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    { -5,  0, 3 } } );

      try {
         const ST sym( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-symmetric SymmetricMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (SymmetricMatrix)
   {
      test_ = "Row-major SymmetricMatrix conversion constructor (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > sym1;
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      const ST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   // Default constructor (StaticMatrix)
   {
      test_ = "Column-major SymmetricMatrix default constructor (StaticMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > sym;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 0UL );
   }

   // Default constructor (HybridMatrix)
   {
      test_ = "Column-major SymmetricMatrix default constructor (HybridMatrix)";

      blaze::SymmetricMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::columnMajor> > sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Default constructor (DynamicMatrix)
   {
      test_ = "Column-major SymmetricMatrix default constructor (DynamicMatrix)";

      OST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   // Size constructor (HybridMatrix)
   {
      test_ = "Column-major SymmetricMatrix size constructor (HybridMatrix)";

      blaze::SymmetricMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::columnMajor> > sym( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkNonZeros( sym, 0UL );
   }

   // Size constructor (DynamicMatrix)
   {
      test_ = "Column-major SymmetricMatrix size constructor (DynamicMatrix)";

      OST sym( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Column-major list initialization
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Column-major SymmetricMatrix initializer list constructor (complete list)";

      const OST sym{ { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Column-major SymmetricMatrix initializer list constructor (incomplete list)";

      const OST sym{ { 1, 2, 3 }, { 2, 4 }, { 3, 0, 6 } };

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major array initialization
   //=====================================================================================

   // Dynamic array initialization constructor
   {
      test_ = "Column-major SymmetricMatrix dynamic array initialization constructor";

      std::unique_ptr<int[]> array( new int[9] );
      array[0] = 1;
      array[1] = 2;
      array[2] = 3;
      array[3] = 2;
      array[4] = 4;
      array[5] = 0;
      array[6] = 3;
      array[7] = 0;
      array[8] = 6;
      const OST sym( 3UL, array.get() );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Static array initialization constructor
   {
      test_ = "Column-major SymmetricMatrix static array initialization constructor";

      const int array[3][3] = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };
      const OST sym( array );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major custom matrix constructors
   //=====================================================================================

   // Custom matrix constructor (ElementType*, size_t)
   {
      test_ = "Column-major SymmetricMatrix custom matrix constructor (ElementType*, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[5UL] );
      memory[1] = 1;
      memory[2] = 2;
      memory[3] = 2;
      memory[4] = 1;
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( memory.get()+1UL, 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 ||
          sym(1,0) != 2 || sym(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, size_t)
   {
      test_ = "Column-major SymmetricMatrix custom matrix constructor (ElementType*, size_t, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[11UL] );
      memory[1] = 1;
      memory[2] = 2;
      memory[6] = 2;
      memory[7] = 1;
      const blaze::SymmetricMatrix<UnalignedUnpadded> sym( memory.get()+1UL, 2UL, 5UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkCapacity( sym, 4UL );
      checkNonZeros( sym, 4UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 ||
          sym(1,0) != 2 || sym(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Column-major SymmetricMatrix copy constructor (0x0)";

      const OST sym1;
      const OST sym2( sym1 );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Column-major SymmetricMatrix copy constructor (3x3)";

      OST sym1( 3UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      const OST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Column-major SymmetricMatrix move constructor (0x0)";

      OST sym1;
      OST sym2( std::move( sym1 ) );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Column-major SymmetricMatrix move constructor (3x3)";

      OST sym1( 3UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      OST sym2( std::move( sym1 ) );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Column-major SymmetricMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat;
      const OST sym( mat );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Conversion constructor (symmetric)
   {
      test_ = "Column-major SymmetricMatrix conversion constructor (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       {  7,  0, 3 } } );

      const OST sym( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 0 ||
          sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-symmetric)
   {
      test_ = "Column-major SymmetricMatrix conversion constructor (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       { -5,  0, 3 } } );

      try {
         const OST sym( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-symmetric SymmetricMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (SymmetricMatrix)
   {
      test_ = "Column-major SymmetricMatrix conversion constructor (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > sym1;
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      const OST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the SymmetricMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseScalarTest::testAssignment()
{
   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Row-major SymmetricMatrix initializer list assignment";

      ST sym;
      sym = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Row-major SymmetricMatrix initializer list assignment";

      ST sym;
      sym = { { 1, 2, 3 }, { 2, 4 }, { 3, 0, 6 } };

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array assignment
   //=====================================================================================

   // Array assignment
   {
      test_ = "Row-major SymmetricMatrix array assignment";

      const int array[3][3] = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };
      ST sym;
      sym = array;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Row-major SymmetricMatrix copy assignment (0x0)";

      ST sym1, sym2;

      sym2 = sym1;

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Row-major SymmetricMatrix copy assignment (3x3)";

      ST sym1( 3UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Row-major SymmetricMatrix move assignment (0x0)";

      ST sym1, sym2;

      sym2 = std::move( sym1 );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Row-major SymmetricMatrix move assignment (3x3)";

      ST sym1( 3UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      ST sym2;
      sym2 = std::move( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major SymmetricMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      ST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Row-major/row-major dense matrix assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix assignment (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    {  7,  0, 3 } } );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 0 ||
          sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix assignment (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       {  7,  0, 3 } } );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 0 ||
          sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix assignment (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    { -5,  0, 3 } } );

      try {
         ST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix assignment (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       { -5,  0, 3 } } );

      try {
         ST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > sym1;
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > sym1;
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major SymmetricMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat;

      ST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Row-major/row-major sparse matrix assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
      mat.insert( 1UL, 2UL, 0 );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 0 ||
          sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
      mat.insert( 1UL, 2UL, 0 );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 0 ||
          sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) = -5;
      mat(2,2) =  3;

      try {
         ST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) = -5;
      mat(2,2) =  3;

      try {
         ST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > sym1( 3UL, 7UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > sym1( 3UL, 7UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Column-major SymmetricMatrix initializer list assignment";

      OST sym;
      sym = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Column-major SymmetricMatrix initializer list assignment";

      OST sym;
      sym = { { 1, 2, 3 }, { 2, 4 }, { 3, 0, 6 } };

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major array assignment
   //=====================================================================================

   // Array assignment
   {
      test_ = "Column-major SymmetricMatrix array assignment";

      const int array[3][3] = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };
      OST sym;
      sym = array;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) != 1 || sym(0,1) != 2 || sym(0,2) != 3 ||
          sym(1,0) != 2 || sym(1,1) != 4 || sym(1,2) != 0 ||
          sym(2,0) != 3 || sym(2,1) != 0 || sym(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Column-major SymmetricMatrix copy assignment (0x0)";

      OST sym1, sym2;

      sym2 = sym1;

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Column-major SymmetricMatrix copy assignment (3x3)";

      OST sym1( 3UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Column-major SymmetricMatrix move assignment (0x0)";

      OST sym1, sym2;

      sym2 = std::move( sym1 );

      checkRows    ( sym2, 0UL );
      checkColumns ( sym2, 0UL );
      checkNonZeros( sym2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Column-major SymmetricMatrix move assignment (3x3)";

      OST sym1( 3UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      OST sym2;
      sym2 = std::move( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major SymmetricMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat;

      OST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Column-major/row-major dense matrix assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix assignment (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    {  7,  0, 3 } } );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 0 ||
          sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix assignment (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       {  7,  0, 3 } } );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 0 ||
          sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix assignment (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    { -5,  0, 3 } } );

      try {
         OST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix assignment (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       { -5,  0, 3 } } );

      try {
         OST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > sym1;
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > sym1;
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major SymmetricMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat;

      OST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Column-major/row-major sparse matrix assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
      mat.insert( 1UL, 2UL, 0 );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 0 ||
          sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
      mat.insert( 1UL, 2UL, 0 );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 7 ||
          sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 0 ||
          sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) = -5;
      mat(2,2) =  3;

      try {
         OST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) = -5;
      mat(2,2) =  3;

      try {
         OST sym;
         sym = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > sym1( 3UL, 7UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > sym1( 3UL, 7UL );
      sym1(0,0) =  1;
      sym1(0,1) = -4;
      sym1(0,2) =  7;
      sym1(1,1) =  2;
      sym1(2,2) =  3;

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -4 || sym2(0,2) != 7 ||
          sym2(1,0) != -4 || sym2(1,1) !=  2 || sym2(1,2) != 0 ||
          sym2(2,0) !=  7 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseScalarTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   // Row-major/row-major dense matrix addition assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix addition assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -6 || sym(0,2) != 13 ||
          sym(1,0) != -6 || sym(1,1) !=  5 || sym(1,2) !=  0 ||
          sym(2,0) != 13 || sym(2,1) !=  0 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix addition assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -6 || sym(0,2) != 13 ||
          sym(1,0) != -6 || sym(1,1) !=  5 || sym(1,2) !=  0 ||
          sym(2,0) != 13 || sym(2,1) !=  0 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix addition assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix addition assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix addition assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix addition assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -6 || sym2(0,2) != 13 ||
          sym2(1,0) != -6 || sym2(1,1) !=  5 || sym2(1,2) !=  0 ||
          sym2(2,0) != 13 || sym2(2,1) !=  0 || sym2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix addition assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -6 || sym2(0,2) != 13 ||
          sym2(1,0) != -6 || sym2(1,1) !=  5 || sym2(1,2) !=  0 ||
          sym2(2,0) != 13 || sym2(2,1) !=  0 || sym2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix addition assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix addition assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -6 || sym(0,2) != 13 ||
          sym(1,0) != -6 || sym(1,1) !=  5 || sym(1,2) !=  0 ||
          sym(2,0) != 13 || sym(2,1) !=  0 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix addition assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -6 || sym(0,2) != 13 ||
          sym(1,0) != -6 || sym(1,1) !=  5 || sym(1,2) !=  0 ||
          sym(2,0) != 13 || sym(2,1) !=  0 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix addition assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix addition assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix addition assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix addition assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix addition assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix addition assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -6 || sym2(0,2) != 13 ||
          sym2(1,0) != -6 || sym2(1,1) !=  5 || sym2(1,2) !=  0 ||
          sym2(2,0) != 13 || sym2(2,1) !=  0 || sym2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix addition assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -6 || sym2(0,2) != 13 ||
          sym2(1,0) != -6 || sym2(1,1) !=  5 || sym2(1,2) !=  0 ||
          sym2(2,0) != 13 || sym2(2,1) !=  0 || sym2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   // Column-major/row-major dense matrix addition assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix addition assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -6 || sym(0,2) != 13 ||
          sym(1,0) != -6 || sym(1,1) !=  5 || sym(1,2) !=  0 ||
          sym(2,0) != 13 || sym(2,1) !=  0 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix addition assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -6 || sym(0,2) != 13 ||
          sym(1,0) != -6 || sym(1,1) !=  5 || sym(1,2) !=  0 ||
          sym(2,0) != 13 || sym(2,1) !=  0 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix addition assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix addition assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix addition assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix addition assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -6 || sym2(0,2) != 13 ||
          sym2(1,0) != -6 || sym2(1,1) !=  5 || sym2(1,2) !=  0 ||
          sym2(2,0) != 13 || sym2(2,1) !=  0 || sym2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix addition assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -6 || sym2(0,2) != 13 ||
          sym2(1,0) != -6 || sym2(1,1) !=  5 || sym2(1,2) !=  0 ||
          sym2(2,0) != 13 || sym2(2,1) !=  0 || sym2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix addition assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix addition assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -6 || sym(0,2) != 13 ||
          sym(1,0) != -6 || sym(1,1) !=  5 || sym(1,2) !=  0 ||
          sym(2,0) != 13 || sym(2,1) !=  0 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix addition assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -6 || sym(0,2) != 13 ||
          sym(1,0) != -6 || sym(1,1) !=  5 || sym(1,2) !=  0 ||
          sym(2,0) != 13 || sym(2,1) !=  0 || sym(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix addition assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix addition assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix addition assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix addition assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix addition assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix addition assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -6 || sym2(0,2) != 13 ||
          sym2(1,0) != -6 || sym2(1,1) !=  5 || sym2(1,2) !=  0 ||
          sym2(2,0) != 13 || sym2(2,1) !=  0 || sym2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix addition assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -6 || sym2(0,2) != 13 ||
          sym2(1,0) != -6 || sym2(1,1) !=  5 || sym2(1,2) !=  0 ||
          sym2(2,0) != 13 || sym2(2,1) !=  0 || sym2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseScalarTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major dense matrix subtraction assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -2 || sym(0,2) != 1 ||
          sym(1,0) != -2 || sym(1,1) != -1 || sym(1,2) != 0 ||
          sym(2,0) !=  1 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -2 || sym(0,2) != 1 ||
          sym(1,0) != -2 || sym(1,1) != -1 || sym(1,2) != 0 ||
          sym(2,0) !=  1 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix subtraction assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix subtraction assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -2 || sym2(0,2) != 1 ||
          sym2(1,0) != -2 || sym2(1,1) != -1 || sym2(1,2) != 0 ||
          sym2(2,0) !=  1 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix subtraction assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -2 || sym2(0,2) != 1 ||
          sym2(1,0) != -2 || sym2(1,1) != -1 || sym2(1,2) != 0 ||
          sym2(2,0) !=  1 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix subtraction assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix subtraction assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -2 || sym(0,2) != 1 ||
          sym(1,0) != -2 || sym(1,1) != -1 || sym(1,2) != 0 ||
          sym(2,0) !=  1 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix subtraction assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -2 || sym(0,2) != 1 ||
          sym(1,0) != -2 || sym(1,1) != -1 || sym(1,2) != 0 ||
          sym(2,0) !=  1 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix subtraction assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix subtraction assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix subtraction assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix subtraction assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -2 || sym2(0,2) != 1 ||
          sym2(1,0) != -2 || sym2(1,1) != -1 || sym2(1,2) != 0 ||
          sym2(2,0) !=  1 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix subtraction assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -2 || sym2(0,2) != 1 ||
          sym2(1,0) != -2 || sym2(1,1) != -1 || sym2(1,2) != 0 ||
          sym2(2,0) !=  1 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major dense matrix subtraction assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -2 || sym(0,2) != 1 ||
          sym(1,0) != -2 || sym(1,1) != -1 || sym(1,2) != 0 ||
          sym(2,0) !=  1 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -2 || sym(0,2) != 1 ||
          sym(1,0) != -2 || sym(1,1) != -1 || sym(1,2) != 0 ||
          sym(2,0) !=  1 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix subtraction assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix subtraction assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -2 || sym2(0,2) != 1 ||
          sym2(1,0) != -2 || sym2(1,1) != -1 || sym2(1,2) != 0 ||
          sym2(2,0) !=  1 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix subtraction assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -2 || sym2(0,2) != 1 ||
          sym2(1,0) != -2 || sym2(1,1) != -1 || sym2(1,2) != 0 ||
          sym2(2,0) !=  1 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix subtraction assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix subtraction assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -2 || sym(0,2) != 1 ||
          sym(1,0) != -2 || sym(1,1) != -1 || sym(1,2) != 0 ||
          sym(2,0) !=  1 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix subtraction assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  1 || sym(0,1) != -2 || sym(0,2) != 1 ||
          sym(1,0) != -2 || sym(1,1) != -1 || sym(1,2) != 0 ||
          sym(2,0) !=  1 || sym(2,1) !=  0 || sym(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix subtraction assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix subtraction assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix subtraction assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix subtraction assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -2 || sym2(0,2) != 1 ||
          sym2(1,0) != -2 || sym2(1,1) != -1 || sym2(1,2) != 0 ||
          sym2(2,0) !=  1 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix subtraction assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  1 || sym2(0,1) != -2 || sym2(0,2) != 1 ||
          sym2(1,0) != -2 || sym2(1,1) != -1 || sym2(1,2) != 0 ||
          sym2(2,0) !=  1 || sym2(2,1) !=  0 || sym2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix Schur product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment operators of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseScalarTest::testSchurAssign()
{
   //=====================================================================================
   // Row-major dense matrix Schur product assignment
   //=====================================================================================

   // Row-major/row-major dense matrix Schur product assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( sym(0,0) !=  0 || sym(0,1) != 8 || sym(0,2) != 42 ||
          sym(1,0) !=  8 || sym(1,1) != 6 || sym(1,2) !=  0 ||
          sym(2,0) != 42 || sym(2,1) != 0 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix Schur product assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( sym(0,0) !=  0 || sym(0,1) != 8 || sym(0,2) != 42 ||
          sym(1,0) !=  8 || sym(1,1) != 6 || sym(1,2) !=  0 ||
          sym(2,0) != 42 || sym(2,1) != 0 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix Schur product assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix Schur product assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( sym2(0,0) !=  0 || sym2(0,1) != 8 || sym2(0,2) != 42 ||
          sym2(1,0) !=  8 || sym2(1,1) != 6 || sym2(1,2) !=  0 ||
          sym2(2,0) != 42 || sym2(2,1) != 0 || sym2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix Schur product assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix Schur product assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( sym2(0,0) !=  0 || sym2(0,1) != 8 || sym2(0,2) != 42 ||
          sym2(1,0) !=  8 || sym2(1,1) != 6 || sym2(1,2) !=  0 ||
          sym2(2,0) != 42 || sym2(2,1) != 0 || sym2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix Schur product assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix Schur product assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix Schur product assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( sym(0,0) !=  0 || sym(0,1) != 8 || sym(0,2) != 42 ||
          sym(1,0) !=  8 || sym(1,1) != 6 || sym(1,2) !=  0 ||
          sym(2,0) != 42 || sym(2,1) != 0 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix Schur product assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix Schur product assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( sym(0,0) !=  0 || sym(0,1) != 8 || sym(0,2) != 42 ||
          sym(1,0) !=  8 || sym(1,1) != 6 || sym(1,2) !=  0 ||
          sym(2,0) != 42 || sym(2,1) != 0 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix Schur product assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix Schur product assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix Schur product assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix Schur product assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( sym2(0,0) !=  0 || sym2(0,1) != 8 || sym2(0,2) != 42 ||
          sym2(1,0) !=  8 || sym2(1,1) != 6 || sym2(1,2) !=  0 ||
          sym2(2,0) != 42 || sym2(2,1) != 0 || sym2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix Schur product assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix Schur product assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( sym2(0,0) !=  0 || sym2(0,1) != 8 || sym2(0,2) != 42 ||
          sym2(1,0) !=  8 || sym2(1,1) != 6 || sym2(1,2) !=  0 ||
          sym2(2,0) != 42 || sym2(2,1) != 0 || sym2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix Schur product assignment
   //=====================================================================================

   // Column-major/row-major dense matrix Schur product assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( sym(0,0) !=  0 || sym(0,1) != 8 || sym(0,2) != 42 ||
          sym(1,0) !=  8 || sym(1,1) != 6 || sym(1,2) !=  0 ||
          sym(2,0) != 42 || sym(2,1) != 0 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix Schur product assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( sym(0,0) !=  0 || sym(0,1) != 8 || sym(0,2) != 42 ||
          sym(1,0) !=  8 || sym(1,1) != 6 || sym(1,2) !=  0 ||
          sym(2,0) != 42 || sym(2,1) != 0 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix Schur product assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix Schur product assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( sym2(0,0) !=  0 || sym2(0,1) != 8 || sym2(0,2) != 42 ||
          sym2(1,0) !=  8 || sym2(1,1) != 6 || sym2(1,2) !=  0 ||
          sym2(2,0) != 42 || sym2(2,1) != 0 || sym2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix Schur product assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix Schur product assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( sym2(0,0) !=  0 || sym2(0,1) != 8 || sym2(0,2) != 42 ||
          sym2(1,0) !=  8 || sym2(1,1) != 6 || sym2(1,2) !=  0 ||
          sym2(2,0) != 42 || sym2(2,1) != 0 || sym2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix Schur product assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix Schur product assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix Schur product assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( sym(0,0) !=  0 || sym(0,1) != 8 || sym(0,2) != 42 ||
          sym(1,0) !=  8 || sym(1,1) != 6 || sym(1,2) !=  0 ||
          sym(2,0) != 42 || sym(2,1) != 0 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix Schur product assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix Schur product assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( sym(0,0) !=  0 || sym(0,1) != 8 || sym(0,2) != 42 ||
          sym(1,0) !=  8 || sym(1,1) != 6 || sym(1,2) !=  0 ||
          sym(2,0) != 42 || sym(2,1) != 0 || sym(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix Schur product assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix Schur product assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix Schur product assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix Schur product assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( sym2(0,0) !=  0 || sym2(0,1) != 8 || sym2(0,2) != 42 ||
          sym2(1,0) !=  8 || sym2(1,1) != 6 || sym2(1,2) !=  0 ||
          sym2(2,0) != 42 || sym2(2,1) != 0 || sym2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix Schur product assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix Schur product assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > sym1( 3UL, 5UL );
      sym1(0,1) = -2;
      sym1(0,2) =  6;
      sym1(1,1) =  3;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( sym2(0,0) !=  0 || sym2(0,1) != 8 || sym2(0,2) != 42 ||
          sym2(1,0) !=  8 || sym2(1,1) != 6 || sym2(1,2) !=  0 ||
          sym2(2,0) != 42 || sym2(2,1) != 0 || sym2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the SymmetricMatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the SymmetricMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseScalarTest::testMultAssign()
{
   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major dense matrix multiplication assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  2 || sym(0,1) != -8 || sym(0,2) != 14 ||
          sym(1,0) != -8 || sym(1,1) !=  4 || sym(1,2) !=  0 ||
          sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  2 || sym(0,1) != -8 || sym(0,2) != 14 ||
          sym(1,0) != -8 || sym(1,1) !=  4 || sym(1,2) !=  0 ||
          sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix multiplication assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix multiplication assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix multiplication assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix multiplication assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,0) = 2;
      sym1(1,1) = 2;
      sym1(2,2) = 2;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  2 || sym2(0,1) != -8 || sym2(0,2) != 14 ||
          sym2(1,0) != -8 || sym2(1,1) !=  4 || sym2(1,2) !=  0 ||
          sym2(2,0) != 14 || sym2(2,1) !=  0 || sym2(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix multiplication assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix multiplication assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,0) = 2;
      sym1(1,1) = 2;
      sym1(2,2) = 2;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  2 || sym2(0,1) != -8 || sym2(0,2) != 14 ||
          sym2(1,0) != -8 || sym2(1,1) !=  4 || sym2(1,2) !=  0 ||
          sym2(2,0) != 14 || sym2(2,1) !=  0 || sym2(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix multiplication assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix multiplication assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  2 || sym(0,1) != -8 || sym(0,2) != 14 ||
          sym(1,0) != -8 || sym(1,1) !=  4 || sym(1,2) !=  0 ||
          sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix multiplication assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  2 || sym(0,1) != -8 || sym(0,2) != 14 ||
          sym(1,0) != -8 || sym(1,1) !=  4 || sym(1,2) !=  0 ||
          sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix multiplication assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix multiplication assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix multiplication assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix multiplication assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      ST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix multiplication assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix multiplication assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > sym1( 3UL, 3UL );
      sym1(0,0) = 2;
      sym1(1,1) = 2;
      sym1(2,2) = 2;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  2 || sym2(0,1) != -8 || sym2(0,2) != 14 ||
          sym2(1,0) != -8 || sym2(1,1) !=  4 || sym2(1,2) !=  0 ||
          sym2(2,0) != 14 || sym2(2,1) !=  0 || sym2(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix multiplication assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix multiplication assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > sym1( 3UL, 3UL );
      sym1(0,0) = 2;
      sym1(1,1) = 2;
      sym1(2,2) = 2;

      ST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  2 || sym2(0,1) != -8 || sym2(0,2) != 14 ||
          sym2(1,0) != -8 || sym2(1,1) !=  4 || sym2(1,2) !=  0 ||
          sym2(2,0) != 14 || sym2(2,1) !=  0 || sym2(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major dense matrix multiplication assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  2 || sym(0,1) != -8 || sym(0,2) != 14 ||
          sym(1,0) != -8 || sym(1,1) !=  4 || sym(1,2) !=  0 ||
          sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  2 || sym(0,1) != -8 || sym(0,2) != 14 ||
          sym(1,0) != -8 || sym(1,1) !=  4 || sym(1,2) !=  0 ||
          sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix multiplication assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix multiplication assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix multiplication assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix multiplication assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,0) = 2;
      sym1(1,1) = 2;
      sym1(2,2) = 2;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  2 || sym2(0,1) != -8 || sym2(0,2) != 14 ||
          sym2(1,0) != -8 || sym2(1,1) !=  4 || sym2(1,2) !=  0 ||
          sym2(2,0) != 14 || sym2(2,1) !=  0 || sym2(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix multiplication assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix multiplication assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,0) = 2;
      sym1(1,1) = 2;
      sym1(2,2) = 2;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  2 || sym2(0,1) != -8 || sym2(0,2) != 14 ||
          sym2(1,0) != -8 || sym2(1,1) !=  4 || sym2(1,2) !=  0 ||
          sym2(2,0) != 14 || sym2(2,1) !=  0 || sym2(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix multiplication assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix multiplication assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  2 || sym(0,1) != -8 || sym(0,2) != 14 ||
          sym(1,0) != -8 || sym(1,1) !=  4 || sym(1,2) !=  0 ||
          sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix multiplication assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,0) = 2;
      mat(1,1) = 2;
      mat(2,2) = 2;
      mat.insert( 1UL, 2UL, 0 );

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 9UL );
      checkNonZeros( sym, 7UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 2UL );

      if( sym(0,0) !=  2 || sym(0,1) != -8 || sym(0,2) != 14 ||
          sym(1,0) != -8 || sym(1,1) !=  4 || sym(1,2) !=  0 ||
          sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix multiplication assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix multiplication assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix multiplication assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix multiplication assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OST sym( 3UL );
      sym(0,0) =  1;
      sym(0,1) = -4;
      sym(0,2) =  7;
      sym(1,1) =  2;
      sym(2,2) =  3;

      try {
         sym *= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix multiplication assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix multiplication assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > sym1( 3UL, 3UL );
      sym1(0,0) = 2;
      sym1(1,1) = 2;
      sym1(2,2) = 2;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  2 || sym2(0,1) != -8 || sym2(0,2) != 14 ||
          sym2(1,0) != -8 || sym2(1,1) !=  4 || sym2(1,2) !=  0 ||
          sym2(2,0) != 14 || sym2(2,1) !=  0 || sym2(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix multiplication assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix multiplication assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > sym1( 3UL, 3UL );
      sym1(0,0) = 2;
      sym1(1,1) = 2;
      sym1(2,2) = 2;

      OST sym2( 3UL );
      sym2(0,0) =  1;
      sym2(0,1) = -4;
      sym2(0,2) =  7;
      sym2(1,1) =  2;
      sym2(2,2) =  3;

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 9UL );
      checkNonZeros( sym2, 7UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 2UL );

      if( sym2(0,0) !=  2 || sym2(0,1) != -8 || sym2(0,2) != 14 ||
          sym2(1,0) != -8 || sym2(1,1) !=  4 || sym2(1,2) !=  0 ||
          sym2(2,0) != 14 || sym2(2,1) !=  0 || sym2(2,2) !=  6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n(  2 -8 14 )\n( -8  4  0 )\n( 14  0  6 )\n";
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
   std::cout << "   Running SymmetricMatrix dense scalar test (part 1)..." << std::endl;

   try
   {
      RUN_SYMMETRICMATRIX_DENSESCALAR_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SymmetricMatrix dense scalar test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
