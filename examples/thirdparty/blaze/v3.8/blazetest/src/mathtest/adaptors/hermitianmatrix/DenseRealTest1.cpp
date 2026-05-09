//=================================================================================================
/*!
//  \file src/mathtest/adaptors/hermitianmatrix/DenseRealTest1.cpp
//  \brief Source file for the HermitianMatrix dense real test (part 1)
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
#include <blazetest/mathtest/adaptors/hermitianmatrix/DenseRealTest.h>

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
/*!\brief Constructor for the HermitianMatrix dense test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseRealTest::DenseRealTest()
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testSchurAssign();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the HermitianMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the HermitianMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseRealTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor (StaticMatrix)
   {
      test_ = "Row-major HermitianMatrix default constructor (StaticMatrix)";

      const blaze::HermitianMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > herm;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 0UL );
   }

   // Default constructor (HybridMatrix)
   {
      test_ = "Row-major HermitianMatrix default constructor (HybridMatrix)";

      const blaze::HermitianMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::rowMajor> > herm;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Default constructor (DynamicMatrix)
   {
      test_ = "Row-major HermitianMatrix default constructor (DynamicMatrix)";

      const HT herm;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   // Size constructor (HybridMatrix)
   {
      test_ = "Row-major HermitianMatrix size constructor (HybridMatrix)";

      const blaze::HermitianMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::rowMajor> > herm( 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkCapacity( herm, 4UL );
      checkNonZeros( herm, 0UL );
   }

   // Size constructor (DynamicMatrix)
   {
      test_ = "Row-major HermitianMatrix size constructor (DynamicMatrix)";

      const HT herm( 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkCapacity( herm, 4UL );
      checkNonZeros( herm, 0UL );
   }


   //=====================================================================================
   // Row-major list initialization
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Row-major HermitianMatrix initializer list constructor (complete list)";

      const HT herm{ { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Row-major HermitianMatrix initializer list constructor (incomplete list)";

      const HT herm{ { 1, 2, 3 }, { 2, 4 }, { 3, 0, 6 } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array initialization
   //=====================================================================================

   // Dynamic array initialization constructor
   {
      test_ = "Row-major HermitianMatrix dynamic array initialization constructor";

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
      const HT herm( 3UL, array.get() );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Static array initialization constructor
   {
      test_ = "Row-major HermitianMatrix static array initialization constructor";

      const int array[3][3] = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };
      const HT herm( array );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major custom matrix constructors
   //=====================================================================================

   // Custom matrix constructor (ElementType*, size_t)
   {
      test_ = "Row-major HermitianMatrix custom matrix constructor (ElementType*, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[5UL] );
      memory[1] = 1;
      memory[2] = 2;
      memory[3] = 2;
      memory[4] = 1;
      const blaze::HermitianMatrix<UnalignedUnpadded> herm( memory.get()+1UL, 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkCapacity( herm, 4UL );
      checkNonZeros( herm, 4UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 ||
          herm(1,0) != 2 || herm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, size_t)
   {
      test_ = "Row-major HermitianMatrix custom matrix constructor (ElementType*, size_t, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[11UL] );
      memory[1] = 1;
      memory[2] = 2;
      memory[6] = 2;
      memory[7] = 1;
      const blaze::HermitianMatrix<UnalignedUnpadded> herm( memory.get()+1UL, 2UL, 5UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkCapacity( herm, 4UL );
      checkNonZeros( herm, 4UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 ||
          herm(1,0) != 2 || herm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Row-major HermitianMatrix copy constructor (0x0)";

      const HT herm1;
      const HT herm2( herm1 );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Row-major HermitianMatrix copy constructor (3x3)";

      HT herm1( 3UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      const HT herm2( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Row-major HermitianMatrix move constructor (0x0)";

      HT herm1;
      HT herm2( std::move( herm1 ) );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Row-major HermitianMatrix move constructor (3x3)";

      HT herm1( 3UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      HT herm2( std::move( herm1 ) );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Row-major HermitianMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;
      const HT herm( mat );

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Conversion constructor (symmetric)
   {
      test_ = "Row-major HermitianMatrix conversion constructor (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    {  7,  0, 3 } } );

      const HT herm( mat );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) !=  1 || herm(0,1) != -4 || herm(0,2) != 7 ||
          herm(1,0) != -4 || herm(1,1) !=  2 || herm(1,2) != 0 ||
          herm(2,0) !=  7 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-symmetric)
   {
      test_ = "Row-major HermitianMatrix conversion constructor (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    { -5,  0, 3 } } );

      try {
         const HT herm( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-symmetric HermitianMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (HermitianMatrix)
   {
      test_ = "Row-major HermitianMatrix conversion constructor (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > herm1;
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      const HT herm2( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   // Default constructor (StaticMatrix)
   {
      test_ = "Column-major HermitianMatrix default constructor (StaticMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > herm;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 0UL );
   }

   // Default constructor (HybridMatrix)
   {
      test_ = "Column-major HermitianMatrix default constructor (HybridMatrix)";

      blaze::HermitianMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::columnMajor> > herm;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Default constructor (DynamicMatrix)
   {
      test_ = "Column-major HermitianMatrix default constructor (DynamicMatrix)";

      OHT herm;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   // Size constructor (HybridMatrix)
   {
      test_ = "Column-major HermitianMatrix size constructor (HybridMatrix)";

      blaze::HermitianMatrix< blaze::HybridMatrix<int,3UL,3UL,blaze::columnMajor> > herm( 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkNonZeros( herm, 0UL );
   }

   // Size constructor (DynamicMatrix)
   {
      test_ = "Column-major HermitianMatrix size constructor (DynamicMatrix)";

      OHT herm( 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkNonZeros( herm, 0UL );
   }


   //=====================================================================================
   // Column-major list initialization
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Column-major HermitianMatrix initializer list constructor (complete list)";

      const OHT herm{ { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Column-major HermitianMatrix initializer list constructor (incomplete list)";

      const OHT herm{ { 1, 2, 3 }, { 2, 4 }, { 3, 0, 6 } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major array initialization
   //=====================================================================================

   // Dynamic array initialization constructor
   {
      test_ = "Column-major HermitianMatrix dynamic array initialization constructor";

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
      const OHT herm( 3UL, array.get() );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Static array initialization constructor
   {
      test_ = "Column-major HermitianMatrix static array initialization constructor";

      const int array[3][3] = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };
      const OHT herm( array );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major custom matrix constructors
   //=====================================================================================

   // Custom matrix constructor (ElementType*, size_t)
   {
      test_ = "Column-major HermitianMatrix custom matrix constructor (ElementType*, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[5UL] );
      memory[1] = 1;
      memory[2] = 2;
      memory[3] = 2;
      memory[4] = 1;
      const blaze::HermitianMatrix<UnalignedUnpadded> herm( memory.get()+1UL, 2UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkCapacity( herm, 4UL );
      checkNonZeros( herm, 4UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 ||
          herm(1,0) != 2 || herm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Custom matrix constructor (ElementType*, size_t, size_t)
   {
      test_ = "Column-major HermitianMatrix custom matrix constructor (ElementType*, size_t, size_t)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[11UL] );
      memory[1] = 1;
      memory[2] = 2;
      memory[6] = 2;
      memory[7] = 1;
      const blaze::HermitianMatrix<UnalignedUnpadded> herm( memory.get()+1UL, 2UL, 5UL );

      checkRows    ( herm, 2UL );
      checkColumns ( herm, 2UL );
      checkCapacity( herm, 4UL );
      checkNonZeros( herm, 4UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 ||
          herm(1,0) != 2 || herm(1,1) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 )\n( 2 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   // Copy constructor (0x0)
   {
      test_ = "Column-major HermitianMatrix copy constructor (0x0)";

      const OHT herm1;
      const OHT herm2( herm1 );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Copy constructor (3x3)
   {
      test_ = "Column-major HermitianMatrix copy constructor (3x3)";

      OHT herm1( 3UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      const OHT herm2( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   // Move constructor (0x0)
   {
      test_ = "Column-major HermitianMatrix move constructor (0x0)";

      OHT herm1;
      OHT herm2( std::move( herm1 ) );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Move constructor (3x3)
   {
      test_ = "Column-major HermitianMatrix move constructor (3x3)";

      OHT herm1( 3UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      OHT herm2( std::move( herm1 ) );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Column-major HermitianMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat;
      const OHT herm( mat );

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Conversion constructor (symmetric)
   {
      test_ = "Column-major HermitianMatrix conversion constructor (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       {  7,  0, 3 } } );

      const OHT herm( mat );

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) !=  1 || herm(0,1) != -4 || herm(0,2) != 7 ||
          herm(1,0) != -4 || herm(1,1) !=  2 || herm(1,2) != 0 ||
          herm(2,0) !=  7 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-symmetric)
   {
      test_ = "Column-major HermitianMatrix conversion constructor (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       { -5,  0, 3 } } );

      try {
         const OHT herm( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of non-symmetric HermitianMatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Conversion constructor (HermitianMatrix)
   {
      test_ = "Column-major HermitianMatrix conversion constructor (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > herm1;
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      const OHT herm2( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the HermitianMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the HermitianMatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseRealTest::testAssignment()
{
   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Row-major HermitianMatrix initializer list assignment (complete list)";

      HT herm;
      herm = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Row-major HermitianMatrix initializer list assignment (incomplete list)";

      HT herm;
      herm = { { 1, 2, 3 }, { 2, 4 }, { 3, 0, 6 } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array assignment
   //=====================================================================================

   // Array assignment
   {
      test_ = "Row-major HermitianMatrix array assignment";

      const int array[3][3] = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };
      HT herm;
      herm = array;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Row-major HermitianMatrix copy assignment (0x0)";

      HT herm1, herm2;

      herm2 = herm1;

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Row-major HermitianMatrix copy assignment (3x3)";

      HT herm1( 3UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      HT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Row-major HermitianMatrix move assignment (0x0)";

      HT herm1, herm2;

      herm2 = std::move( herm1 );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Row-major HermitianMatrix move assignment (3x3)";

      HT herm1( 3UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      HT herm2;
      herm2 = std::move( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major HermitianMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      HT herm;
      herm = mat;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Row-major/row-major dense matrix assignment (symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix assignment (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    {  7,  0, 3 } } );

      HT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) !=  1 || herm(0,1) != -4 || herm(0,2) != 7 ||
          herm(1,0) != -4 || herm(1,1) !=  2 || herm(1,2) != 0 ||
          herm(2,0) !=  7 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix assignment (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       {  7,  0, 3 } } );

      HT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) !=  1 || herm(0,1) != -4 || herm(0,2) != 7 ||
          herm(1,0) != -4 || herm(1,1) !=  2 || herm(1,2) != 0 ||
          herm(2,0) !=  7 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix assignment (non-symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix assignment (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    { -5,  0, 3 } } );

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix assignment (non-symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix assignment (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       { -5,  0, 3 } } );

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > herm1;
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      HT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > herm1;
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      HT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major HermitianMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat;

      HT herm;
      herm = mat;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Row-major/row-major sparse matrix assignment (symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
      mat.insert( 1UL, 2UL, 0 );

      HT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) !=  1 || herm(0,1) != -4 || herm(0,2) != 7 ||
          herm(1,0) != -4 || herm(1,1) !=  2 || herm(1,2) != 0 ||
          herm(2,0) !=  7 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
      mat.insert( 1UL, 2UL, 0 );

      HT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) !=  1 || herm(0,1) != -4 || herm(0,2) != 7 ||
          herm(1,0) != -4 || herm(1,1) !=  2 || herm(1,2) != 0 ||
          herm(2,0) !=  7 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) = -5;
      mat(2,2) =  3;

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) = -5;
      mat(2,2) =  3;

      try {
         HT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > herm1( 3UL, 7UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      HT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > herm1( 3UL, 7UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      HT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   // Complete initializer list
   {
      test_ = "Column-major HermitianMatrix initializer list assignment (complete list)";

      OHT herm;
      herm = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Incomplete initializer list
   {
      test_ = "Column-major HermitianMatrix initializer list assignment (incomplete list)";

      OHT herm;
      herm = { { 1, 2, 3 }, { 2, 4 }, { 3, 0, 6 } };

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major array assignment
   //=====================================================================================

   // Array assignment
   {
      test_ = "Column-major HermitianMatrix array assignment";

      const int array[3][3] = { { 1, 2, 3 }, { 2, 4, 0 }, { 3, 0, 6 } };
      OHT herm;
      herm = array;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) != 1 || herm(0,1) != 2 || herm(0,2) != 3 ||
          herm(1,0) != 2 || herm(1,1) != 4 || herm(1,2) != 0 ||
          herm(2,0) != 3 || herm(2,1) != 0 || herm(2,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 2 4 0 )\n( 3 0 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   // Copy assignment (0x0)
   {
      test_ = "Column-major HermitianMatrix copy assignment (0x0)";

      OHT herm1, herm2;

      herm2 = herm1;

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Copy assignment (3x3)
   {
      test_ = "Column-major HermitianMatrix copy assignment (3x3)";

      OHT herm1( 3UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      OHT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   // Move assignment (0x0)
   {
      test_ = "Column-major HermitianMatrix move assignment (0x0)";

      OHT herm1, herm2;

      herm2 = std::move( herm1 );

      checkRows    ( herm2, 0UL );
      checkColumns ( herm2, 0UL );
      checkNonZeros( herm2, 0UL );
   }

   // Move assignment (3x3)
   {
      test_ = "Column-major HermitianMatrix move assignment (3x3)";

      OHT herm1( 3UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      OHT herm2;
      herm2 = std::move( herm1 );

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major HermitianMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat;

      OHT herm;
      herm = mat;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Column-major/row-major dense matrix assignment (symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix assignment (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    {  7,  0, 3 } } );

      OHT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) !=  1 || herm(0,1) != -4 || herm(0,2) != 7 ||
          herm(1,0) != -4 || herm(1,1) !=  2 || herm(1,2) != 0 ||
          herm(2,0) !=  7 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix assignment (symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       {  7,  0, 3 } } );

      OHT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) !=  1 || herm(0,1) != -4 || herm(0,2) != 7 ||
          herm(1,0) != -4 || herm(1,1) !=  2 || herm(1,2) != 0 ||
          herm(2,0) !=  7 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix assignment (non-symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix assignment (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> mat( { {  1, -4, 7 },
                                                                    { -4,  2, 0 },
                                                                    { -5,  0, 3 } } );

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix assignment (non-symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix assignment (non-symmetric)";

      const blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> mat( { {  1, -4, 7 },
                                                                       { -4,  2, 0 },
                                                                       { -5,  0, 3 } } );

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::rowMajor> > herm1;
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      OHT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::StaticMatrix<int,3UL,3UL,blaze::columnMajor> > herm1;
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      OHT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major HermitianMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat;

      OHT herm;
      herm = mat;

      checkRows    ( herm, 0UL );
      checkColumns ( herm, 0UL );
      checkNonZeros( herm, 0UL );
   }

   // Column-major/row-major sparse matrix assignment (symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
      mat.insert( 1UL, 2UL, 0 );

      OHT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) !=  1 || herm(0,1) != -4 || herm(0,2) != 7 ||
          herm(1,0) != -4 || herm(1,1) !=  2 || herm(1,2) != 0 ||
          herm(2,0) !=  7 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 8UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) =  7;
      mat(2,2) =  3;
      mat.insert( 1UL, 2UL, 0 );

      OHT herm;
      herm = mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( herm(0,0) !=  1 || herm(0,1) != -4 || herm(0,2) != 7 ||
          herm(1,0) != -4 || herm(1,1) !=  2 || herm(1,2) != 0 ||
          herm(2,0) !=  7 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) = -5;
      mat(2,2) =  3;

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) =  1;
      mat(0,1) = -4;
      mat(0,2) =  7;
      mat(1,0) = -4;
      mat(1,1) =  2;
      mat(2,0) = -5;
      mat(2,2) =  3;

      try {
         OHT herm;
         herm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > herm1( 3UL, 7UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      OHT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > herm1( 3UL, 7UL );
      herm1(0,0) =  1;
      herm1(0,1) = -4;
      herm1(0,2) =  7;
      herm1(1,1) =  2;
      herm1(2,2) =  3;

      OHT herm2;
      herm2 = herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkNonZeros( herm2, 7UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -4 || herm2(0,2) != 7 ||
          herm2(1,0) != -4 || herm2(1,1) !=  2 || herm2(1,2) != 0 ||
          herm2(2,0) !=  7 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -4  7 )\n( -4  2  0 )\n(  7  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the HermitianMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseRealTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   // Row-major/row-major dense matrix addition assignment (symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix addition assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -6 || herm(0,2) != 13 ||
          herm(1,0) != -6 || herm(1,1) !=  5 || herm(1,2) !=  0 ||
          herm(2,0) != 13 || herm(2,1) !=  0 || herm(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix addition assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -6 || herm(0,2) != 13 ||
          herm(1,0) != -6 || herm(1,1) !=  5 || herm(1,2) !=  0 ||
          herm(2,0) != 13 || herm(2,1) !=  0 || herm(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix addition assignment (non-symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix addition assignment (non-symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      HT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -6 || herm2(0,2) != 13 ||
          herm2(1,0) != -6 || herm2(1,1) !=  5 || herm2(1,2) !=  0 ||
          herm2(2,0) != 13 || herm2(2,1) !=  0 || herm2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      OHT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -6 || herm2(0,2) != 13 ||
          herm2(1,0) != -6 || herm2(1,1) !=  5 || herm2(1,2) !=  0 ||
          herm2(2,0) != 13 || herm2(2,1) !=  0 || herm2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix addition assignment (symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix addition assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -6 || herm(0,2) != 13 ||
          herm(1,0) != -6 || herm(1,1) !=  5 || herm(1,2) !=  0 ||
          herm(2,0) != 13 || herm(2,1) !=  0 || herm(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix addition assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -6 || herm(0,2) != 13 ||
          herm(1,0) != -6 || herm(1,1) !=  5 || herm(1,2) !=  0 ||
          herm(2,0) != 13 || herm(2,1) !=  0 || herm(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix addition assignment (non-symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix addition assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix addition assignment (non-symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix addition assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix addition assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix addition assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -6 || herm2(0,2) != 13 ||
          herm2(1,0) != -6 || herm2(1,1) !=  5 || herm2(1,2) !=  0 ||
          herm2(2,0) != 13 || herm2(2,1) !=  0 || herm2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix addition assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix addition assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -6 || herm2(0,2) != 13 ||
          herm2(1,0) != -6 || herm2(1,1) !=  5 || herm2(1,2) !=  0 ||
          herm2(2,0) != 13 || herm2(2,1) !=  0 || herm2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   // Column-major/row-major dense matrix addition assignment (symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix addition assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -6 || herm(0,2) != 13 ||
          herm(1,0) != -6 || herm(1,1) !=  5 || herm(1,2) !=  0 ||
          herm(2,0) != 13 || herm(2,1) !=  0 || herm(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix addition assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -6 || herm(0,2) != 13 ||
          herm(1,0) != -6 || herm(1,1) !=  5 || herm(1,2) !=  0 ||
          herm(2,0) != 13 || herm(2,1) !=  0 || herm(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix addition assignment (non-symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix addition assignment (non-symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      HT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -6 || herm2(0,2) != 13 ||
          herm2(1,0) != -6 || herm2(1,1) !=  5 || herm2(1,2) !=  0 ||
          herm2(2,0) != 13 || herm2(2,1) !=  0 || herm2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix addition assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix addition assignment (HermitianMatrix)";

      OHT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -6 || herm2(0,2) != 13 ||
          herm2(1,0) != -6 || herm2(1,1) !=  5 || herm2(1,2) !=  0 ||
          herm2(2,0) != 13 || herm2(2,1) !=  0 || herm2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix addition assignment (symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix addition assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -6 || herm(0,2) != 13 ||
          herm(1,0) != -6 || herm(1,1) !=  5 || herm(1,2) !=  0 ||
          herm(2,0) != 13 || herm(2,1) !=  0 || herm(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix addition assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm += mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -6 || herm(0,2) != 13 ||
          herm(1,0) != -6 || herm(1,1) !=  5 || herm(1,2) !=  0 ||
          herm(2,0) != 13 || herm(2,1) !=  0 || herm(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix addition assignment (non-symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix addition assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix addition assignment (non-symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix addition assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix addition assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix addition assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -6 || herm2(0,2) != 13 ||
          herm2(1,0) != -6 || herm2(1,1) !=  5 || herm2(1,2) !=  0 ||
          herm2(2,0) != 13 || herm2(2,1) !=  0 || herm2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix addition assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix addition assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 += herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -6 || herm2(0,2) != 13 ||
          herm2(1,0) != -6 || herm2(1,1) !=  5 || herm2(1,2) !=  0 ||
          herm2(2,0) != 13 || herm2(2,1) !=  0 || herm2(2,2) !=  3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -6 13 )\n( -6  5  0 )\n( 13  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the HermitianMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseRealTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major dense matrix subtraction assignment (symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -2 || herm(0,2) != 1 ||
          herm(1,0) != -2 || herm(1,1) != -1 || herm(1,2) != 0 ||
          herm(2,0) !=  1 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -2 || herm(0,2) != 1 ||
          herm(1,0) != -2 || herm(1,1) != -1 || herm(1,2) != 0 ||
          herm(2,0) !=  1 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      HT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -2 || herm2(0,2) != 1 ||
          herm2(1,0) != -2 || herm2(1,1) != -1 || herm2(1,2) != 0 ||
          herm2(2,0) !=  1 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      OHT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -2 || herm2(0,2) != 1 ||
          herm2(1,0) != -2 || herm2(1,1) != -1 || herm2(1,2) != 0 ||
          herm2(2,0) !=  1 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix subtraction assignment (symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix subtraction assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -2 || herm(0,2) != 1 ||
          herm(1,0) != -2 || herm(1,1) != -1 || herm(1,2) != 0 ||
          herm(2,0) !=  1 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix subtraction assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -2 || herm(0,2) != 1 ||
          herm(1,0) != -2 || herm(1,1) != -1 || herm(1,2) != 0 ||
          herm(2,0) !=  1 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix subtraction assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix subtraction assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix subtraction assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -2 || herm2(0,2) != 1 ||
          herm2(1,0) != -2 || herm2(1,1) != -1 || herm2(1,2) != 0 ||
          herm2(2,0) !=  1 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix subtraction assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -2 || herm2(0,2) != 1 ||
          herm2(1,0) != -2 || herm2(1,1) != -1 || herm2(1,2) != 0 ||
          herm2(2,0) !=  1 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major dense matrix subtraction assignment (symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -2 || herm(0,2) != 1 ||
          herm(1,0) != -2 || herm(1,1) != -1 || herm(1,2) != 0 ||
          herm(2,0) !=  1 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -2 || herm(0,2) != 1 ||
          herm(1,0) != -2 || herm(1,1) != -1 || herm(1,2) != 0 ||
          herm(2,0) !=  1 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      HT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -2 || herm2(0,2) != 1 ||
          herm2(1,0) != -2 || herm2(1,1) != -1 || herm2(1,2) != 0 ||
          herm2(2,0) !=  1 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix subtraction assignment (HermitianMatrix)";

      OHT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -2 || herm2(0,2) != 1 ||
          herm2(1,0) != -2 || herm2(1,1) != -1 || herm2(1,2) != 0 ||
          herm2(2,0) !=  1 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix subtraction assignment (symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix subtraction assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -2 || herm(0,2) != 1 ||
          herm(1,0) != -2 || herm(1,1) != -1 || herm(1,2) != 0 ||
          herm(2,0) !=  1 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix subtraction assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm -= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 7UL );
      checkNonZeros( herm, 0UL, 3UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 2UL );

      if( herm(0,0) !=  1 || herm(0,1) != -2 || herm(0,2) != 1 ||
          herm(1,0) != -2 || herm(1,1) != -1 || herm(1,2) != 0 ||
          herm(2,0) !=  1 || herm(2,1) !=  0 || herm(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix subtraction assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix subtraction assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix subtraction assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -2 || herm2(0,2) != 1 ||
          herm2(1,0) != -2 || herm2(1,1) != -1 || herm2(1,2) != 0 ||
          herm2(2,0) !=  1 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix subtraction assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix subtraction assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 -= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 7UL );
      checkNonZeros( herm2, 0UL, 3UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 2UL );

      if( herm2(0,0) !=  1 || herm2(0,1) != -2 || herm2(0,2) != 1 ||
          herm2(1,0) != -2 || herm2(1,1) != -1 || herm2(1,2) != 0 ||
          herm2(2,0) !=  1 || herm2(2,1) !=  0 || herm2(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1 -2  1 )\n( -2 -1  0 )\n(  1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the HermitianMatrix Schur product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment operators of the HermitianMatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseRealTest::testSchurAssign()
{
   //=====================================================================================
   // Row-major dense matrix Schur product assignment
   //=====================================================================================

   // Row-major/row-major dense matrix Schur product assignment (symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) !=  0 || herm(0,1) != 8 || herm(0,2) != 42 ||
          herm(1,0) !=  8 || herm(1,1) != 6 || herm(1,2) !=  0 ||
          herm(2,0) != 42 || herm(2,1) != 0 || herm(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix Schur product assignment (symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) !=  0 || herm(0,1) != 8 || herm(0,2) != 42 ||
          herm(1,0) !=  8 || herm(1,1) != 6 || herm(1,2) !=  0 ||
          herm(2,0) != 42 || herm(2,1) != 0 || herm(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major dense matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix dense matrix Schur product assignment (HermitianMatrix)";

      HT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) !=  0 || herm2(0,1) != 8 || herm2(0,2) != 42 ||
          herm2(1,0) !=  8 || herm2(1,1) != 6 || herm2(1,2) !=  0 ||
          herm2(2,0) != 42 || herm2(2,1) != 0 || herm2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix dense matrix Schur product assignment (HermitianMatrix)";

      OHT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) !=  0 || herm2(0,1) != 8 || herm2(0,2) != 42 ||
          herm2(1,0) !=  8 || herm2(1,1) != 6 || herm2(1,2) !=  0 ||
          herm2(2,0) != 42 || herm2(2,1) != 0 || herm2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix Schur product assignment
   //=====================================================================================

   // Row-major/row-major sparse matrix Schur product assignment (symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix Schur product assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) !=  0 || herm(0,1) != 8 || herm(0,2) != 42 ||
          herm(1,0) !=  8 || herm(1,1) != 6 || herm(1,2) !=  0 ||
          herm(2,0) != 42 || herm(2,1) != 0 || herm(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix Schur product assignment (symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix Schur product assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) !=  0 || herm(0,1) != 8 || herm(0,2) != 42 ||
          herm(1,0) !=  8 || herm(1,1) != 6 || herm(1,2) !=  0 ||
          herm(2,0) != 42 || herm(2,1) != 0 || herm(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix Schur product assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major sparse matrix Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix Schur product assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      HT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/row-major sparse matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Row-major/row-major HermitianMatrix sparse matrix Schur product assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) !=  0 || herm2(0,1) != 8 || herm2(0,2) != 42 ||
          herm2(1,0) !=  8 || herm2(1,1) != 6 || herm2(1,2) !=  0 ||
          herm2(2,0) != 42 || herm2(2,1) != 0 || herm2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Row-major/column-major HermitianMatrix sparse matrix Schur product assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      HT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) !=  0 || herm2(0,1) != 8 || herm2(0,2) != 42 ||
          herm2(1,0) !=  8 || herm2(1,1) != 6 || herm2(1,2) !=  0 ||
          herm2(2,0) != 42 || herm2(2,1) != 0 || herm2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix Schur product assignment
   //=====================================================================================

   // Column-major/row-major dense matrix Schur product assignment (symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) !=  0 || herm(0,1) != 8 || herm(0,2) != 42 ||
          herm(1,0) !=  8 || herm(1,1) != 6 || herm(1,2) !=  0 ||
          herm(2,0) != 42 || herm(2,1) != 0 || herm(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix Schur product assignment (symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) !=  0 || herm(0,1) != 8 || herm(0,2) != 42 ||
          herm(1,0) !=  8 || herm(1,1) != 6 || herm(1,2) !=  0 ||
          herm(2,0) != 42 || herm(2,1) != 0 || herm(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 0 );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major dense matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix dense matrix Schur product assignment (HermitianMatrix)";

      HT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) !=  0 || herm2(0,1) != 8 || herm2(0,2) != 42 ||
          herm2(1,0) !=  8 || herm2(1,1) != 6 || herm2(1,2) !=  0 ||
          herm2(2,0) != 42 || herm2(2,1) != 0 || herm2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix dense matrix Schur product assignment (HermitianMatrix)";

      OHT herm1( 3UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) !=  0 || herm2(0,1) != 8 || herm2(0,2) != 42 ||
          herm2(1,0) !=  8 || herm2(1,1) != 6 || herm2(1,2) !=  0 ||
          herm2(2,0) != 42 || herm2(2,1) != 0 || herm2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix Schur product assignment
   //=====================================================================================

   // Column-major/row-major sparse matrix Schur product assignment (symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix Schur product assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) !=  0 || herm(0,1) != 8 || herm(0,2) != 42 ||
          herm(1,0) !=  8 || herm(1,1) != 6 || herm(1,2) !=  0 ||
          herm(2,0) != 42 || herm(2,1) != 0 || herm(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix Schur product assignment (symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix Schur product assignment (symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 6UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,0) = -2;
      mat(1,1) =  3;
      mat(2,0) =  6;
      mat.insert( 1UL, 2UL, 0 );

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      herm %= mat;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkCapacity( herm, 9UL );
      checkNonZeros( herm, 5UL );
      checkNonZeros( herm, 0UL, 2UL );
      checkNonZeros( herm, 1UL, 2UL );
      checkNonZeros( herm, 2UL, 1UL );

      if( herm(0,0) !=  0 || herm(0,1) != 8 || herm(0,2) != 42 ||
          herm(1,0) !=  8 || herm(1,1) != 6 || herm(1,2) !=  0 ||
          herm(2,0) != 42 || herm(2,1) != 0 || herm(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix Schur product assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major sparse matrix Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix Schur product assignment (non-symmetric)";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 3UL, 3UL, 4UL );
      mat(0,1) = -2;
      mat(0,2) =  6;
      mat(1,1) =  3;
      mat(2,0) =  6;

      OHT herm( 3UL );
      herm(0,0) =  1;
      herm(0,1) = -4;
      herm(0,2) =  7;
      herm(1,1) =  2;
      herm(2,2) =  3;

      try {
         herm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/row-major sparse matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Column-major/row-major HermitianMatrix sparse matrix Schur product assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) !=  0 || herm2(0,1) != 8 || herm2(0,2) != 42 ||
          herm2(1,0) !=  8 || herm2(1,1) != 6 || herm2(1,2) !=  0 ||
          herm2(2,0) != 42 || herm2(2,1) != 0 || herm2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix Schur product assignment (HermitianMatrix)
   {
      test_ = "Column-major/column-major HermitianMatrix sparse matrix Schur product assignment (HermitianMatrix)";

      blaze::HermitianMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > herm1( 3UL, 5UL );
      herm1(0,1) = -2;
      herm1(0,2) =  6;
      herm1(1,1) =  3;

      OHT herm2( 3UL );
      herm2(0,0) =  1;
      herm2(0,1) = -4;
      herm2(0,2) =  7;
      herm2(1,1) =  2;
      herm2(2,2) =  3;

      herm2 %= herm1;

      checkRows    ( herm2, 3UL );
      checkColumns ( herm2, 3UL );
      checkCapacity( herm2, 9UL );
      checkNonZeros( herm2, 5UL );
      checkNonZeros( herm2, 0UL, 2UL );
      checkNonZeros( herm2, 1UL, 2UL );
      checkNonZeros( herm2, 2UL, 1UL );

      if( herm2(0,0) !=  0 || herm2(0,1) != 8 || herm2(0,2) != 42 ||
          herm2(1,0) !=  8 || herm2(1,1) != 6 || herm2(1,2) !=  0 ||
          herm2(2,0) != 42 || herm2(2,1) != 0 || herm2(2,2) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << herm2 << "\n"
             << "   Expected result:\n(  1  8 42 )\n(  8  6  0 )\n( 42  0  0 )\n";
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
   std::cout << "   Running HermitianMatrix dense real test (part 1)..." << std::endl;

   try
   {
      RUN_HERMITIANMATRIX_DENSEREAL_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during HermitianMatrix dense real test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
