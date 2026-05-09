//=================================================================================================
/*!
//  \file src/mathtest/adaptors/symmetricmatrix/SparseNonScalarTest1.cpp
//  \brief Source file for the SymmetricMatrix sparse non-scalar test (part 1)
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blazetest/mathtest/adaptors/symmetricmatrix/SparseNonScalarTest.h>

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
/*!\brief Constructor for the SymmetricMatrix sparse non-scalar test.
//
// \exception std::runtime_error Operation error detected.
*/
SparseNonScalarTest::SparseNonScalarTest()
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
void SparseNonScalarTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor (CompressedMatrix)
   {
      test_ = "Row-major SymmetricMatrix default constructor (CompressedMatrix)";

      const ST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   // Size constructor (CompressedMatrix)
   {
      test_ = "Row-major SymmetricMatrix size constructor (CompressedMatrix)";

      const ST sym( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkNonZeros( sym, 0UL );
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
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      const ST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 7UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) ( 0 ) )\n"
                                     "( (  7 ) (  0 ) ( 3 ) )\n";
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
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2( std::move( sym1 ) );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 7UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) ( 0 ) )\n"
                                     "( (  7 ) (  0 ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Row-major SymmetricMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<VT,blaze::rowMajor> mat;
      const ST sym( mat );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Conversion constructor (symmetric)
   {
      test_ = "Row-major SymmetricMatrix conversion constructor (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      const ST sym( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 7UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-symmetric)
   {
      test_ = "Row-major SymmetricMatrix conversion constructor (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

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

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      const ST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 7UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   // Default constructor (CompressedMatrix)
   {
      test_ = "Column-major SymmetricMatrix default constructor (CompressedMatrix)";

      const OST sym;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   // Size constructor (CompressedMatrix)
   {
      test_ = "Column-major SymmetricMatrix size constructor (CompressedMatrix)";

      const OST sym( 2UL );

      checkRows    ( sym, 2UL );
      checkColumns ( sym, 2UL );
      checkNonZeros( sym, 0UL );
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
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      const OST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 7UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
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
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2( std::move( sym1 ) );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 7UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major conversion constructor
   //=====================================================================================

   // Conversion constructor (0x0)
   {
      test_ = "Column-major SymmetricMatrix conversion constructor (0x0)";

      const blaze::DynamicMatrix<VT,blaze::columnMajor> mat;
      const OST sym( mat );

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Conversion constructor (symmetric)
   {
      test_ = "Column-major SymmetricMatrix conversion constructor (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      const OST sym( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkCapacity( sym, 7UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Conversion constructor (non-symmetric)
   {
      test_ = "Column-major SymmetricMatrix conversion constructor (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

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

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      const OST sym2( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkCapacity( sym2, 7UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
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
void SparseNonScalarTest::testAssignment()
{
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
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
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
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = std::move( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major SymmetricMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<VT,blaze::rowMajor> mat;

      ST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Row-major/row-major dense matrix assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

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

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

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

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Row-major SymmetricMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<VT,blaze::rowMajor> mat;

      ST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Row-major/row-major sparse matrix assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

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

      blaze::CompressedMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

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

      blaze::SymmetricMatrix< blaze::CompressedMatrix<VT,blaze::rowMajor> > sym1( 3UL, 5UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major sparse matrix assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix sparse matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<VT,blaze::columnMajor> > sym1( 3UL, 5UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      ST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix computation assignment
   //=====================================================================================

   // Row-major/row-major dense matrix computation assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix computation assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major dense matrix computation assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix computation assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      ST sym;
      sym = eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major dense matrix computation assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix dense matrix computation assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         ST sym;
         sym = eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major dense matrix computation assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix dense matrix computation assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         ST sym;
         sym = eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
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
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
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
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = std::move( sym1 );

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major SymmetricMatrix dense matrix assignment (0x0)";

      const blaze::DynamicMatrix<VT,blaze::columnMajor> mat;

      OST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Column-major/row-major dense matrix assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

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

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

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

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> > sym1;
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   // Conversion assignment (0x0)
   {
      test_ = "Column-major SymmetricMatrix sparse matrix assignment (0x0)";

      const blaze::CompressedMatrix<VT,blaze::columnMajor> mat;

      OST sym;
      sym = mat;

      checkRows    ( sym, 0UL );
      checkColumns ( sym, 0UL );
      checkNonZeros( sym, 0UL );
   }

   // Column-major/row-major sparse matrix assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix assignment (symmetric)";

      blaze::CompressedMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major sparse matrix assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix sparse matrix assignment (non-symmetric)";

      blaze::CompressedMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

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

      blaze::CompressedMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, 7UL );
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

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

      blaze::SymmetricMatrix< blaze::CompressedMatrix<VT,blaze::rowMajor> > sym1( 3UL, 5UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major sparse matrix assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix sparse matrix assignment (SymmetricMatrix)";

      blaze::SymmetricMatrix< blaze::CompressedMatrix<VT,blaze::columnMajor> > sym1( 3UL, 5UL );
      sym1(0,0) = vec(  1 );
      sym1(0,1) = vec( -4 );
      sym1(0,2) = vec(  7 );
      sym1(1,1) = vec(  2 );
      sym1(2,2) = vec(  3 );

      OST sym2;
      sym2 = sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 7UL );

      if( sym2(0,0) != vec(  1 ) || sym2(0,1) != vec( -4 )  || sym2(0,2) != vec( 7 )   ||
          sym2(1,0) != vec( -4 ) || sym2(1,1) != vec(  2 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  7 ) || !isDefault( sym2(2,1) ) || sym2(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major computation assignment
   //=====================================================================================

   // Column-major/row-major dense matrix computation assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix computation assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major dense matrix computation assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix computation assignment (symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec(  7 );
      mat(2,2) = vec(  3 );

      OST sym;
      sym = eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 7UL );

      if( sym(0,0) != vec(  1 ) || sym(0,1) != vec( -4 )  || sym(0,2) != vec( 7 )   ||
          sym(1,0) != vec( -4 ) || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  7 ) || !isDefault( sym(2,1) ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  1 ) ( -4 ) ( 7 ) )\n"
                                     "( ( -4 ) (  2 ) (   ) )\n"
                                     "( (  7 ) (    ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major dense matrix computation assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix dense matrix computation assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::rowMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         OST sym;
         sym = eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major dense matrix computation assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix dense matrix computation assignment (non-symmetric)";

      blaze::StaticMatrix<VT,3UL,3UL,blaze::columnMajor> mat;
      mat(0,0) = vec(  1 );
      mat(0,1) = vec( -4 );
      mat(0,2) = vec(  7 );
      mat(1,0) = vec( -4 );
      mat(1,1) = vec(  2 );
      mat(2,0) = vec( -5 );
      mat(2,2) = vec(  3 );

      try {
         OST sym;
         sym = eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
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
void SparseNonScalarTest::testAddAssign()
{
   //=====================================================================================
   // Row-major addition assignment
   //=====================================================================================

   // Row-major/row-major addition assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major addition assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major addition assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Row-major/column-major addition assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Row-major/row-major addition assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix addition assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      ST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major addition assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix addition assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      ST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major computation addition assignment
   //=====================================================================================

   // Row-major/row-major computation addition assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix computation addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major computation addition assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix computation addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major computation addition assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix computation addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major computation addition assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix computation addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major addition assignment
   //=====================================================================================

   // Column-major/row-major addition assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major addition assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major addition assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Column-major/column-major addition assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Column-major/row-major addition assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix addition assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      OST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major addition assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix addition assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      OST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 += sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major computation addition assignment
   //=====================================================================================

   // Column-major/row-major computation addition assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix computation addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major computation addition assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix computation addition assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym += eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major computation addition assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix computation addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major computation addition assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix computation addition assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym += eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
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
void SparseNonScalarTest::testSubAssign()
{
   //=====================================================================================
   // Row-major subtraction assignment
   //=====================================================================================

   // Row-major/row-major subtraction assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major subtraction assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Row-major/column-major subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Row-major/row-major subtraction assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix subtraction assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = vec(  2 );
      sym1(0,2) = vec( -6 );
      sym1(1,1) = vec( -3 );

      ST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major subtraction assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix subtraction assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = vec(  2 );
      sym1(0,2) = vec( -6 );
      sym1(1,1) = vec( -3 );

      ST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major computation computation subtraction assignment
   //=====================================================================================

   // Row-major/row-major computation subtraction assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix computation subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major computation subtraction assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix computation subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major computation subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix computation subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major computation subtraction assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix computation subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major subtraction assignment
   //=====================================================================================

   // Column-major/row-major subtraction assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major subtraction assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n( ( -1 ) (  3 ) (   ) )\n( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Column-major/column-major subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Column-major/row-major subtraction assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix subtraction assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = vec(  2 );
      sym1(0,2) = vec( -6 );
      sym1(1,1) = vec( -3 );

      OST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major subtraction assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix subtraction assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = vec(  2 );
      sym1(0,2) = vec( -6 );
      sym1(1,1) = vec( -3 );

      OST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 -= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -1 )  || sym2(0,2) != vec( 8 )   ||
          sym2(1,0) != vec( -1 )  || sym2(1,1) != vec(  3 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec(  8 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major computation subtraction assignment
   //=====================================================================================

   // Column-major/row-major computation subtraction assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix computation subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major computation subtraction assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix computation subtraction assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,0) = vec(  2 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym -= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -1 )  || sym(0,2) != vec( 8 )   ||
          sym(1,0) != vec( -1 )  || sym(1,1) != vec(  3 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec(  8 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -1 ) ( 8 ) )\n"
                                     "( ( -1 ) (  3 ) (   ) )\n"
                                     "( (  8 ) (    ) (   ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major computation subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix computation subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major computation subtraction assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix computation subtraction assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec(  2 );
      mat(0,2) = vec( -6 );
      mat(1,1) = vec( -3 );
      mat(2,0) = vec( -6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym -= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
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
void SparseNonScalarTest::testSchurAssign()
{
   //=====================================================================================
   // Row-major Schur product assignment
   //=====================================================================================

   // Row-major/row-major Schur product assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -2 )  || sym(0,2) != vec( 12 )  ||
          sym(1,0) != vec( -2 )  || sym(1,1) != vec(  0 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 12 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major Schur product assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -2 )  || sym(0,2) != vec( 12 )  ||
          sym(1,0) != vec( -2 )  || sym(1,1) != vec(  0 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 12 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Row-major/column-major Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Row-major/row-major Schur product assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix Schur product assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      ST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -2 )  || sym2(0,2) != vec( 12 )  ||
          sym2(1,0) != vec( -2 )  || sym2(1,1) != vec(  0 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec( 12 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major Schur product assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix Schur product assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      ST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -2 )  || sym2(0,2) != vec( 12 )  ||
          sym2(1,0) != vec( -2 )  || sym2(1,1) != vec(  0 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec( 12 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major computation Schur product assignment
   //=====================================================================================

   // Row-major/row-major computation Schur product assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix computation Schur product assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym %= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -2 )  || sym(0,2) != vec( 12 )  ||
          sym(1,0) != vec( -2 )  || sym(1,1) != vec(  0 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 12 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major computation Schur product assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix computation Schur product assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym %= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -2 )  || sym(0,2) != vec( 12 )  ||
          sym(1,0) != vec( -2 )  || sym(1,1) != vec(  0 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 12 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major computation Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix computation Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym %= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Row-major/column-major computation Schur product assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix computation Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym %= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Column-major Schur product assignment
   //=====================================================================================

   // Column-major/row-major Schur product assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -2 )  || sym(0,2) != vec( 12 )  ||
          sym(1,0) != vec( -2 )  || sym(1,1) != vec(  0 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 12 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major Schur product assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix Schur product assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym %= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -2 )  || sym(0,2) != vec( 12 )  ||
          sym(1,0) != vec( -2 )  || sym(1,1) != vec(  0 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 12 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Column-major/column-major Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

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

   // Column-major/row-major Schur product assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix Schur product assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      OST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -2 )  || sym2(0,2) != vec( 12 )  ||
          sym2(1,0) != vec( -2 )  || sym2(1,1) != vec(  0 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec( 12 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major Schur product assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix Schur product assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,1) = vec( -2 );
      sym1(0,2) = vec(  6 );
      sym1(1,1) = vec(  3 );

      OST sym2( 3UL );
      sym2(0,1) = vec( 1 );
      sym2(0,2) = vec( 2 );
      sym2(1,1) = vec( 0 );

      sym2 %= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 5UL );
      checkNonZeros( sym2, 0UL, 2UL );
      checkNonZeros( sym2, 1UL, 2UL );
      checkNonZeros( sym2, 2UL, 1UL );

      if( !isDefault( sym2(0,0) ) || sym2(0,1) != vec( -2 )  || sym2(0,2) != vec( 12 )  ||
          sym2(1,0) != vec( -2 )  || sym2(1,1) != vec(  0 )  || !isDefault( sym2(1,2) ) ||
          sym2(2,0) != vec( 12 )  || !isDefault( sym2(2,1) ) || !isDefault( sym2(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major computation Schur product assignment
   //=====================================================================================

   // Column-major/row-major computation Schur product assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix computation Schur product assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym %= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -2 )  || sym(0,2) != vec( 12 )  ||
          sym(1,0) != vec( -2 )  || sym(1,1) != vec(  0 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 12 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major computation Schur product assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix computation Schur product assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,0) = vec( -2 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      sym %= eval( mat );

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 5UL );
      checkNonZeros( sym, 0UL, 2UL );
      checkNonZeros( sym, 1UL, 2UL );
      checkNonZeros( sym, 2UL, 1UL );

      if( !isDefault( sym(0,0) ) || sym(0,1) != vec( -2 )  || sym(0,2) != vec( 12 )  ||
          sym(1,0) != vec( -2 )  || sym(1,1) != vec(  0 )  || !isDefault( sym(1,2) ) ||
          sym(2,0) != vec( 12 )  || !isDefault( sym(2,1) ) || !isDefault( sym(2,2) ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (    ) ( -2 ) ( 12 ) )\n"
                                     "( ( -2 ) (  0 ) (    ) )\n"
                                     "( ( 12 ) (    ) (    ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major computation Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix computation Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym %= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric row-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // Column-major/column-major computation Schur product assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix computation Schur product assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,1) = vec( 1 );
      sym(0,2) = vec( 2 );
      sym(1,1) = vec( 0 );

      try {
         sym %= eval( mat );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment of non-symmetric column-major matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
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
void SparseNonScalarTest::testMultAssign()
{
   //=====================================================================================
   // Row-major multiplication assignment
   //=====================================================================================

   // Row-major/row-major multiplication assignment (symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,0) = vec( 2 );
      mat(1,1) = vec( 2 );
      mat(2,2) = vec( 2 );

      ST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec(  2 ) || sym(0,1) != vec( -8 ) || sym(0,2) != vec( 14 ) ||
          sym(1,0) != vec( -8 ) || sym(1,1) != vec(  4 ) || sym(1,2) != vec(  0 ) ||
          sym(2,0) != vec( 14 ) || sym(2,1) != vec(  0 ) || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major multiplication assignment (symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,0) = vec( 2 );
      mat(1,1) = vec( 2 );
      mat(2,2) = vec( 2 );

      ST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec(  2 ) || sym(0,1) != vec( -8 ) || sym(0,2) != vec( 14 ) ||
          sym(1,0) != vec( -8 ) || sym(1,1) != vec(  4 ) || sym(1,2) != vec(  0 ) ||
          sym(2,0) != vec( 14 ) || sym(2,1) != vec(  0 ) || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/row-major multiplication assignment (non-symmetric)
   {
      test_ = "Row-major/row-major SymmetricMatrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

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

   // Row-major/column-major multiplication assignment (non-symmetric)
   {
      test_ = "Row-major/column-major SymmetricMatrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      ST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

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

   // Row-major/row-major multiplication assignment (SymmetricMatrix)
   {
      test_ = "Row-major/row-major SymmetricMatrix multiplication assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,0) = vec( 2 );
      sym1(0,1) = vec( 0 );
      sym1(0,2) = vec( 0 );
      sym1(1,1) = vec( 2 );
      sym1(1,2) = vec( 0 );
      sym1(2,2) = vec( 2 );

      ST sym2( 3UL );
      sym2(0,0) = vec(  1 );
      sym2(0,1) = vec( -4 );
      sym2(0,2) = vec(  7 );
      sym2(1,1) = vec(  2 );
      sym2(1,2) = vec(  0 );
      sym2(2,2) = vec(  3 );

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 9UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 3UL );
      checkNonZeros( sym2, 2UL, 3UL );

      if( sym2(0,0) != vec(  2 ) || sym2(0,1) != vec( -8 ) || sym2(0,2) != vec( 14 ) ||
          sym2(1,0) != vec( -8 ) || sym2(1,1) != vec(  4 ) || sym2(1,2) != vec(  0 ) ||
          sym2(2,0) != vec( 14 ) || sym2(2,1) != vec(  0 ) || sym2(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Row-major/column-major multiplication assignment (SymmetricMatrix)
   {
      test_ = "Row-major/column-major SymmetricMatrix multiplication assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,0) = vec( 2 );
      sym1(0,1) = vec( 0 );
      sym1(0,2) = vec( 0 );
      sym1(1,1) = vec( 2 );
      sym1(1,2) = vec( 0 );
      sym1(2,2) = vec( 2 );

      ST sym2( 3UL );
      sym2(0,0) = vec(  1 );
      sym2(0,1) = vec( -4 );
      sym2(0,2) = vec(  7 );
      sym2(1,1) = vec(  2 );
      sym2(1,2) = vec(  0 );
      sym2(2,2) = vec(  3 );

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 9UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 3UL );
      checkNonZeros( sym2, 2UL, 3UL );

      if( sym2(0,0) != vec(  2 ) || sym2(0,1) != vec( -8 ) || sym2(0,2) != vec( 14 ) ||
          sym2(1,0) != vec( -8 ) || sym2(1,1) != vec(  4 ) || sym2(1,2) != vec(  0 ) ||
          sym2(2,0) != vec( 14 ) || sym2(2,1) != vec(  0 ) || sym2(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major multiplication assignment
   //=====================================================================================

   // Column-major/row-major multiplication assignment (symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,0) = vec( 2 );
      mat(1,1) = vec( 2 );
      mat(2,2) = vec( 2 );

      OST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec(  2 ) || sym(0,1) != vec( -8 ) || sym(0,2) != vec( 14 ) ||
          sym(1,0) != vec( -8 ) || sym(1,1) != vec(  4 ) || sym(1,2) != vec(  0 ) ||
          sym(2,0) != vec( 14 ) || sym(2,1) != vec(  0 ) || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major multiplication assignment (symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix multiplication assignment (symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,0) = vec( 2 );
      mat(1,1) = vec( 2 );
      mat(2,2) = vec( 2 );

      OST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

      sym *= mat;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 9UL );
      checkNonZeros( sym, 0UL, 3UL );
      checkNonZeros( sym, 1UL, 3UL );
      checkNonZeros( sym, 2UL, 3UL );

      if( sym(0,0) != vec(  2 ) || sym(0,1) != vec( -8 ) || sym(0,2) != vec( 14 ) ||
          sym(1,0) != vec( -8 ) || sym(1,1) != vec(  4 ) || sym(1,2) != vec(  0 ) ||
          sym(2,0) != vec( 14 ) || sym(2,1) != vec(  0 ) || sym(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/row-major multiplication assignment (non-symmetric)
   {
      test_ = "Column-major/row-major SymmetricMatrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::rowMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

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

   // Column-major/column-major multiplication assignment (non-symmetric)
   {
      test_ = "Column-major/column-major SymmetricMatrix multiplication assignment (non-symmetric)";

      blaze::DynamicMatrix<VT,blaze::columnMajor> mat( 3UL, 3UL, vec( 0 ) );
      mat(0,1) = vec( -2 );
      mat(0,2) = vec(  6 );
      mat(1,1) = vec(  3 );
      mat(2,0) = vec(  6 );

      OST sym( 3UL );
      sym(0,0) = vec(  1 );
      sym(0,1) = vec( -4 );
      sym(0,2) = vec(  7 );
      sym(1,1) = vec(  2 );
      sym(1,2) = vec(  0 );
      sym(2,2) = vec(  3 );

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

   // Column-major/row-major multiplication assignment (SymmetricMatrix)
   {
      test_ = "Column-major/row-major SymmetricMatrix multiplication assignment (SymmetricMatrix)";

      ST sym1( 3UL );
      sym1(0,0) = vec( 2 );
      sym1(0,1) = vec( 0 );
      sym1(0,2) = vec( 0 );
      sym1(1,1) = vec( 2 );
      sym1(1,2) = vec( 0 );
      sym1(2,2) = vec( 2 );

      OST sym2( 3UL );
      sym2(0,0) = vec(  1 );
      sym2(0,1) = vec( -4 );
      sym2(0,2) = vec(  7 );
      sym2(1,1) = vec(  2 );
      sym2(1,2) = vec(  0 );
      sym2(2,2) = vec(  3 );

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 9UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 3UL );
      checkNonZeros( sym2, 2UL, 3UL );

      if( sym2(0,0) != vec(  2 ) || sym2(0,1) != vec( -8 ) || sym2(0,2) != vec( 14 ) ||
          sym2(1,0) != vec( -8 ) || sym2(1,1) != vec(  4 ) || sym2(1,2) != vec(  0 ) ||
          sym2(2,0) != vec( 14 ) || sym2(2,1) != vec(  0 ) || sym2(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Column-major/column-major multiplication assignment (SymmetricMatrix)
   {
      test_ = "Column-major/column-major SymmetricMatrix multiplication assignment (SymmetricMatrix)";

      OST sym1( 3UL );
      sym1(0,0) = vec( 2 );
      sym1(0,1) = vec( 0 );
      sym1(0,2) = vec( 0 );
      sym1(1,1) = vec( 2 );
      sym1(1,2) = vec( 0 );
      sym1(2,2) = vec( 2 );

      OST sym2( 3UL );
      sym2(0,0) = vec(  1 );
      sym2(0,1) = vec( -4 );
      sym2(0,2) = vec(  7 );
      sym2(1,1) = vec(  2 );
      sym2(1,2) = vec(  0 );
      sym2(2,2) = vec(  3 );

      sym2 *= sym1;

      checkRows    ( sym2, 3UL );
      checkColumns ( sym2, 3UL );
      checkNonZeros( sym2, 9UL );
      checkNonZeros( sym2, 0UL, 3UL );
      checkNonZeros( sym2, 1UL, 3UL );
      checkNonZeros( sym2, 2UL, 3UL );

      if( sym2(0,0) != vec(  2 ) || sym2(0,1) != vec( -8 ) || sym2(0,2) != vec( 14 ) ||
          sym2(1,0) != vec( -8 ) || sym2(1,1) != vec(  4 ) || sym2(1,2) != vec(  0 ) ||
          sym2(2,0) != vec( 14 ) || sym2(2,1) != vec(  0 ) || sym2(2,2) != vec(  6 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sym2 << "\n"
             << "   Expected result:\n( (  2 ) ( -8 ) ( 14 ) )\n"
                                     "( ( -8 ) (  4 ) (  0 ) )\n"
                                     "( ( 14 ) (  0 ) (  6 ) )\n";
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
   std::cout << "   Running SymmetricMatrix sparse non-scalar test (part 1)..." << std::endl;

   try
   {
      RUN_SYMMETRICMATRIX_SPARSENONSCALAR_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SymmetricMatrix sparse non-scalar test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
