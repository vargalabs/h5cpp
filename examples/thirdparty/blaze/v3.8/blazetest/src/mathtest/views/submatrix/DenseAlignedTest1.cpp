//=================================================================================================
/*!
//  \file src/mathtest/views/submatrix/DenseAlignedTest1.cpp
//  \brief Source file for the Submatrix dense aligned test (part 1)
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/CustomMatrix.h>
#include <blaze/util/Memory.h>
#include <blaze/util/policies/Deallocate.h>
#include <blaze/util/typetraits/AlignmentOf.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>
#include <blazetest/mathtest/views/submatrix/DenseAlignedTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace views {

namespace submatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the Submatrix dense aligned test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseAlignedTest::DenseAlignedTest()
   : mat1_ ( 64UL, 64UL )
   , mat2_ ( 64UL, 64UL )
   , tmat1_( 64UL, 64UL )
   , tmat2_( 64UL, 64UL )
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
/*!\brief Test of the Submatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the Submatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testConstructors()
{
   using blaze::submatrix;
   using blaze::aligned;
   using blaze::unaligned;


   //=====================================================================================
   // Row-major submatrix tests
   //=====================================================================================

   {
      test_ = "Row-major Submatrix constructor";

      initialize();

      const size_t alignment = blaze::AlignmentOf<int>::value;

      for( size_t row=0UL; row<mat1_.rows(); row+=alignment ) {
         for( size_t column=0UL; column<mat1_.columns(); column+=alignment ) {
            for( size_t maxm=0UL; ; maxm+=alignment ) {
               for( size_t maxn=0UL; ; maxn+=alignment )
               {
                  const size_t m( blaze::min( maxm, mat1_.rows()-row ) );
                  const size_t n( blaze::min( maxn, mat1_.columns()-column ) );

                  const ASMT sm1 = submatrix<aligned>  ( mat1_, row, column, m, n );
                  const USMT sm2 = submatrix<unaligned>( mat2_, row, column, m, n );

                  if( sm1 != sm2 ) {
                     std::ostringstream oss;
                     oss << " Test: " << test_ << "\n"
                         << " Error: Setup of dense submatrix failed\n"
                         << " Details:\n"
                         << "   Index of first row    = " << row << "\n"
                         << "   Index of first column = " << column << "\n"
                         << "   Number of rows        = " << m << "\n"
                         << "   Number of columns     = " << n << "\n"
                         << "   Submatrix:\n" << sm1 << "\n"
                         << "   Reference:\n" << sm2 << "\n";
                     throw std::runtime_error( oss.str() );
                  }

                  if( column+maxn > mat1_.columns() ) break;
               }

               if( row+maxm > mat1_.rows() ) break;
            }
         }
      }

      try {
         ASMT sm = submatrix<aligned>( mat1_, 0UL, 16UL, 64UL, 49UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         ASMT sm = submatrix<aligned>( mat1_, 16UL, 0UL, 49UL, 64UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         ASMT sm = submatrix<aligned>( mat1_, 80UL, 0UL, 8UL, 8UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         ASMT sm = submatrix<aligned>( mat1_, 0UL, 80UL, 8UL, 8UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      if( blaze::AlignmentOf<int>::value > sizeof(int) )
      {
         try {
            ASMT sm = submatrix<aligned>( mat1_, 8UL, 7UL, 8UL, 8UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of unaligned submatrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Column-major submatrix tests
   //=====================================================================================

   {
      test_ = "Column-major Submatrix constructor";

      initialize();

      const size_t alignment = blaze::AlignmentOf<int>::value;

      for( size_t column=0UL; column<mat1_.columns(); column+=alignment ) {
         for( size_t row=0UL; row<mat1_.rows(); row+=alignment ) {
            for( size_t maxn=0UL; ; maxn+=alignment ) {
               for( size_t maxm=0UL; ; maxm+=alignment )
               {
                  const size_t n( blaze::min( maxn, mat1_.columns()-column ) );
                  const size_t m( blaze::min( maxm, mat1_.rows()-row ) );

                  const AOSMT sm1 = submatrix<aligned>  ( tmat1_, row, column, m, n );
                  const UOSMT sm2 = submatrix<unaligned>( tmat2_, row, column, m, n );

                  if( sm1 != sm2 ) {
                     std::ostringstream oss;
                     oss << " Test: " << test_ << "\n"
                         << " Error: Setup of dense submatrix failed\n"
                         << " Details:\n"
                         << "   Index of first row    = " << row << "\n"
                         << "   Index of first column = " << column << "\n"
                         << "   Number of rows        = " << m << "\n"
                         << "   Number of columns     = " << n << "\n"
                         << "   Submatrix:\n" << sm1 << "\n"
                         << "   Reference:\n" << sm2 << "\n";
                     throw std::runtime_error( oss.str() );
                  }

                  if( row+maxm > mat1_.rows() ) break;
               }

               if( column+maxn > mat1_.columns() ) break;
            }
         }
      }

      try {
         AOSMT sm = submatrix<aligned>( tmat1_, 0UL, 16UL, 64UL, 49UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         AOSMT sm = submatrix<aligned>( tmat1_, 16UL, 0UL, 49UL, 64UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         AOSMT sm = submatrix<aligned>( tmat1_, 80UL, 0UL, 8UL, 8UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      try {
         AOSMT sm = submatrix<aligned>( tmat1_, 0UL, 80UL, 8UL, 8UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds submatrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      if( blaze::AlignmentOf<int>::value > sizeof(int) )
      {
         try {
            AOSMT sm = submatrix<aligned>( tmat1_, 7UL, 8UL, 8UL, 8UL );

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of unaligned submatrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Submatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the Submatrix specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testAssignment()
{
   using blaze::submatrix;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::initializer_list;


   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Row-major Submatrix homogeneous assignment";

      initialize();

      // Assigning to a 8x16 submatrix
      {
         ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
         USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );
         sm1 = 12;
         sm2 = 12;

         checkRows   ( sm1,  8UL );
         checkColumns( sm1, 16UL );
         checkRows   ( sm2,  8UL );
         checkColumns( sm2, 16UL );

         if( sm1 != sm2 || mat1_ != mat2_ ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << sm1 << "\n"
                << "   Expected result:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Assigning to a 16x8 submatrix
      {
         ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 16UL, 8UL );
         USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 16UL, 8UL );
         sm1 = 15;
         sm2 = 15;

         checkRows   ( sm1, 16UL );
         checkColumns( sm1,  8UL );
         checkRows   ( sm2, 16UL );
         checkColumns( sm2,  8UL );

         if( sm1 != sm2 || mat1_ != mat2_ ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << sm1 << "\n"
                << "   Expected result:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   {
      test_ = "Row-major initializer list assignment (complete list)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      initializer_list< initializer_list<int> > list =
         { { 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  13,  14,  15,  16 },
           { 2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24,  26,  28,  30,  32 },
           { 3,  6,  9, 12, 15, 18, 21, 24, 27, 30, 33, 36,  39,  42,  45,  48 },
           { 4,  8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48,  52,  56,  60,  64 },
           { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60,  65,  70,  75,  80 },
           { 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72,  78,  86,  92,  98 },
           { 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84,  91,  98, 105, 112 },
           { 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128 } };

      sm1 = list;
      sm2 = list;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major initializer list assignment (incomplete list)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      initializer_list< initializer_list<int> > list =
         { { 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  13,  14,  15,  16 },
           { 2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24,  26,  28 },
           { 3,  6,  9, 12, 15, 18, 21, 24, 27, 30, 33, 36 },
           { 4,  8, 12, 16, 20, 24, 28, 32, 36, 40 },
           { 5, 10, 15, 20, 25, 30, 35, 40 },
           { 6, 12, 18, 24, 30, 36 },
           { 7, 14, 21, 28 },
           { 8, 16 } };

      sm1 = list;
      sm2 = list;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major Submatrix copy assignment (no aliasing)";

      initialize();

      MT mat1( 64UL, 64UL );
      MT mat2( 64UL, 64UL );
      randomize( mat1, int(randmin), int(randmax) );
      mat2 = mat1;

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );
      sm1 = submatrix<aligned>  ( mat1, 8UL, 16UL, 8UL, 16UL );
      sm2 = submatrix<unaligned>( mat2, 8UL, 16UL, 8UL, 16UL );

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Submatrix copy assignment (aliasing)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );
      sm1 = submatrix<aligned>  ( mat1_, 12UL, 16UL, 8UL, 16UL );
      sm2 = submatrix<unaligned>( mat2_, 12UL, 16UL, 8UL, 16UL );

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix assignment (mixed type)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::DynamicMatrix<short,rowMajor> mat( 8UL, 16UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix assignment (aligned/padded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 8UL, 16UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix assignment (mixed type)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::DynamicMatrix<short,columnMajor> mat( 8UL, 16UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix assignment (aligned/padded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 256UL ) );
      AlignedPadded mat( memory.get(), 8UL, 16UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix assignment";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 8UL, 16UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix assignment";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 8UL, 16UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Column-major Submatrix homogeneous assignment";

      initialize();

      // Assigning to a 8x16 submatrix
      {
         AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 8UL, 16UL );
         UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 8UL, 16UL );
         sm1 = 12;
         sm2 = 12;

         checkRows   ( sm1,  8UL );
         checkColumns( sm1, 16UL );
         checkRows   ( sm2,  8UL );
         checkColumns( sm2, 16UL );

         if( sm1 != sm2 || mat1_ != mat2_ ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << sm1 << "\n"
                << "   Expected result:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Assigning to a 16x8 submatrix
      {
         AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
         UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );
         sm1 = 15;
         sm2 = 15;

         checkRows   ( sm1, 16UL );
         checkColumns( sm1,  8UL );
         checkRows   ( sm2, 16UL );
         checkColumns( sm2,  8UL );

         if( sm1 != sm2 || mat1_ != mat2_ ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << sm1 << "\n"
                << "   Expected result:\n" << sm2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   {
      test_ = "Column-major initializer list assignment (complete list)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      initializer_list< initializer_list<int> > list =
         { {  1,  2,  3,  4,  5,  6,   7,   8 },
           {  2,  4,  6,  8, 10, 12,  14,  16 },
           {  3,  6,  9, 12, 15, 18,  21,  24 },
           {  4,  8, 12, 16, 20, 24,  28,  32 },
           {  5, 10, 15, 20, 25, 30,  35,  40 },
           {  6, 12, 18, 24, 30, 36,  42,  48 },
           {  7, 14, 21, 28, 35, 42,  49,  56 },
           {  8, 16, 24, 32, 40, 48,  56,  64 },
           {  9, 18, 27, 36, 45, 54,  63,  72 },
           { 10, 20, 30, 40, 50, 60,  70,  80 },
           { 11, 22, 33, 44, 55, 66,  77,  88 },
           { 12, 24, 36, 48, 60, 72,  84,  96 },
           { 13, 26, 39, 52, 65, 78,  91, 104 },
           { 14, 28, 42, 56, 70, 84,  98, 112 },
           { 15, 30, 45, 60, 75, 90, 105, 120 },
           { 16, 32, 48, 64, 80, 96, 112, 128 } };

      sm1 = list;
      sm2 = list;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major initializer list assignment (incomplete list)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      initializer_list< initializer_list<int> > list =
         { {  1,  2,  3,  4,  5,  6,   7,   8 },
           {  2,  4,  6,  8, 10, 12,  14 },
           {  3,  6,  9, 12, 15, 18 },
           {  4,  8, 12, 16, 20 },
           {  5, 10, 15, 20 },
           {  6, 12, 18 },
           {  7, 14 },
           {  8 },
           {  9, 18, 27, 36, 45, 54,  63,  72 },
           { 10, 20, 30, 40, 50, 60,  70 },
           { 11, 22, 33, 44, 55, 66 },
           { 12, 24, 36, 48, 60 },
           { 13, 26, 39, 52 },
           { 14, 28, 42 },
           { 15, 30 },
           { 16 } };

      sm1 = list;
      sm2 = list;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major Submatrix copy assignment (no aliasing)";

      initialize();

      OMT mat1( 64UL, 64UL );
      OMT mat2( 64UL, 64UL );
      randomize( mat1, int(randmin), int(randmax) );
      mat2 = mat1;

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );
      sm1 = submatrix<aligned>  ( mat1, 16UL, 8UL, 16UL, 8UL );
      sm2 = submatrix<unaligned>( mat2, 16UL, 8UL, 16UL, 8UL );

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Submatrix copy assignment (aliasing)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );
      sm1 = submatrix<aligned>  ( tmat1_, 16UL, 12UL, 16UL, 8UL );
      sm2 = submatrix<unaligned>( tmat2_, 16UL, 12UL, 16UL, 8UL );

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix assignment (mixed type)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::DynamicMatrix<short,rowMajor> mat( 16UL, 8UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix assignment (aligned/padded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 256UL ) );
      AlignedPadded mat( memory.get(), 16UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 16UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix assignment (mixed type)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::DynamicMatrix<short,columnMajor> mat( 16UL, 8UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix assignment (aligned/padded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 16UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix assignment (unaligned/unpadded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 16UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix assignment";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 16UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix assignment";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 16UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 = mat;
      sm2 = mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Submatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the Submatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testAddAssign()
{
   using blaze::submatrix;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major Submatrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major Submatrix addition assignment (no aliasing)";

      initialize();

      MT mat1( 64UL, 64UL );
      MT mat2( 64UL, 64UL );
      randomize( mat1, int(randmin), int(randmax) );
      mat2 = mat1;

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );
      sm1 += submatrix<aligned>  ( mat1, 8UL, 16UL, 8UL, 16UL );
      sm2 += submatrix<unaligned>( mat2, 8UL, 16UL, 8UL, 16UL );

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Submatrix addition assignment (aliasing)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );
      sm1 += submatrix<aligned>  ( mat1_, 12UL, 16UL, 8UL, 16UL );
      sm2 += submatrix<unaligned>( mat2_, 12UL, 16UL, 8UL, 16UL );

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix addition assignment (mixed type)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::DynamicMatrix<short,rowMajor> mat( 8UL, 16UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix addition assignment (aligned/padded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 8UL, 16UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix addition assignment (mixed type)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::DynamicMatrix<short,columnMajor> mat( 8UL, 16UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix addition assignment (aligned/padded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 256UL ) );
      AlignedPadded mat( memory.get(), 8UL, 16UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix addition assignment";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 8UL, 16UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix addition assignment";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 8UL, 16UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Submatrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major Submatrix addition assignment (no aliasing)";

      initialize();

      OMT mat1( 64UL, 64UL );
      OMT mat2( 64UL, 64UL );
      randomize( mat1, int(randmin), int(randmax) );
      mat2 = mat1;

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );
      sm1 += submatrix<aligned>  ( mat1, 16UL, 8UL, 16UL, 8UL );
      sm2 += submatrix<unaligned>( mat2, 16UL, 8UL, 16UL, 8UL );

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Submatrix addition assignment (aliasing)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );
      sm1 += submatrix<aligned>  ( tmat1_, 16UL, 12UL, 16UL, 8UL );
      sm2 += submatrix<unaligned>( tmat2_, 16UL, 12UL, 16UL, 8UL );

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix addition assignment (mixed type)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::DynamicMatrix<short,rowMajor> mat( 16UL, 8UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix addition assignment (aligned/padded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 256UL ) );
      AlignedPadded mat( memory.get(), 16UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 16UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix addition assignment (mixed type)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::DynamicMatrix<short,columnMajor> mat( 16UL, 8UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix addition assignment (aligned/padded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 16UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix addition assignment (unaligned/unpadded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 16UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix addition assignment";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 16UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix addition assignment";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 16UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 += mat;
      sm2 += mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Submatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the Submatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testSubAssign()
{
   using blaze::submatrix;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major Submatrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major Submatrix subtraction assignment (no aliasing)";

      initialize();

      MT mat1( 64UL, 64UL );
      MT mat2( 64UL, 64UL );
      randomize( mat1, int(randmin), int(randmax) );
      mat2 = mat1;

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );
      sm1 -= submatrix<aligned>  ( mat1, 8UL, 16UL, 8UL, 16UL );
      sm2 -= submatrix<unaligned>( mat2, 8UL, 16UL, 8UL, 16UL );

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Submatrix subtraction assignment (aliasing)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );
      sm1 -= submatrix<aligned>  ( mat1_, 12UL, 16UL, 8UL, 16UL );
      sm2 -= submatrix<unaligned>( mat2_, 12UL, 16UL, 8UL, 16UL );

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix subtraction assignment (mixed type)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::DynamicMatrix<short,rowMajor> mat( 8UL, 16UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 8UL, 16UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix subtraction assignment (mixed type)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::DynamicMatrix<short,columnMajor> mat( 8UL, 16UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 256UL ) );
      AlignedPadded mat( memory.get(), 8UL, 16UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix subtraction assignment";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 8UL, 16UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix subtraction assignment";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 8UL, 16UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Submatrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major Submatrix subtraction assignment (no aliasing)";

      initialize();

      OMT mat1( 64UL, 64UL );
      OMT mat2( 64UL, 64UL );
      randomize( mat1, int(randmin), int(randmax) );
      mat2 = mat1;

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );
      sm1 -= submatrix<aligned>  ( mat1, 16UL, 8UL, 16UL, 8UL );
      sm2 -= submatrix<unaligned>( mat2, 16UL, 8UL, 16UL, 8UL );

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Submatrix subtraction assignment (aliasing)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );
      sm1 -= submatrix<aligned>  ( tmat1_, 16UL, 12UL, 16UL, 8UL );
      sm2 -= submatrix<unaligned>( tmat2_, 16UL, 12UL, 16UL, 8UL );

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix subtraction assignment (mixed type)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::DynamicMatrix<short,rowMajor> mat( 16UL, 8UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 256UL ) );
      AlignedPadded mat( memory.get(), 16UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 16UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix subtraction assignment (mixed type)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::DynamicMatrix<short,columnMajor> mat( 16UL, 8UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix subtraction assignment (aligned/padded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 16UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix subtraction assignment (unaligned/unpadded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 16UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix subtraction assignment";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 16UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix subtraction assignment";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 16UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 -= mat;
      sm2 -= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Submatrix Schur product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment operators of the Submatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testSchurAssign()
{
   using blaze::submatrix;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major Submatrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major Submatrix Schur product assignment (no aliasing)";

      initialize();

      MT mat1( 64UL, 64UL );
      MT mat2( 64UL, 64UL );
      randomize( mat1, int(randmin), int(randmax) );
      mat2 = mat1;

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );
      sm1 %= submatrix<aligned>  ( mat1, 8UL, 16UL, 8UL, 16UL );
      sm2 %= submatrix<unaligned>( mat2, 8UL, 16UL, 8UL, 16UL );

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Submatrix Schur product assignment (aliasing)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );
      sm1 %= submatrix<aligned>  ( mat1_, 12UL, 16UL, 8UL, 16UL );
      sm2 %= submatrix<unaligned>( mat2_, 12UL, 16UL, 8UL, 16UL );

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix Schur product assignment (mixed type)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::DynamicMatrix<short,rowMajor> mat( 8UL, 16UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix Schur product assignment (aligned/padded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 8UL, 16UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix Schur product assignment (unaligned/unpadded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix Schur product assignment (mixed type)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::DynamicMatrix<short,columnMajor> mat( 8UL, 16UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix Schur product assignment (aligned/padded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 256UL ) );
      AlignedPadded mat( memory.get(), 8UL, 16UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix Schur product assignment (unaligned/unpadded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix Schur product assignment";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 8UL, 16UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix Schur product assignment";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 8UL, 16UL, 8UL, 16UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 8UL, 16UL, 8UL, 16UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 8UL, 16UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1,  8UL );
      checkColumns( sm1, 16UL );
      checkRows   ( sm2,  8UL );
      checkColumns( sm2, 16UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Submatrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major Submatrix Schur product assignment (no aliasing)";

      initialize();

      OMT mat1( 64UL, 64UL );
      OMT mat2( 64UL, 64UL );
      randomize( mat1, int(randmin), int(randmax) );
      mat2 = mat1;

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );
      sm1 %= submatrix<aligned>  ( mat1, 16UL, 8UL, 16UL, 8UL );
      sm2 %= submatrix<unaligned>( mat2, 16UL, 8UL, 16UL, 8UL );

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Submatrix Schur product assignment (aliasing)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );
      sm1 %= submatrix<aligned>  ( tmat1_, 16UL, 12UL, 16UL, 8UL );
      sm2 %= submatrix<unaligned>( tmat2_, 16UL, 12UL, 16UL, 8UL );

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix Schur product assignment (mixed type)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::DynamicMatrix<short,rowMajor> mat( 16UL, 8UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix Schur product assignment (aligned/padded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 256UL ) );
      AlignedPadded mat( memory.get(), 16UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix Schur product assignment (unaligned/unpadded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 16UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix Schur product assignment (mixed type)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::DynamicMatrix<short,columnMajor> mat( 16UL, 8UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix Schur product assignment (aligned/padded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 16UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix Schur product assignment (unaligned/unpadded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[129UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 16UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix Schur product assignment";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 16UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix Schur product assignment";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 8UL, 16UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 8UL, 16UL, 8UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 16UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 %= mat;
      sm2 %= mat;

      checkRows   ( sm1, 16UL );
      checkColumns( sm1,  8UL );
      checkRows   ( sm2, 16UL );
      checkColumns( sm2,  8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Submatrix multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the Submatrix
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testMultAssign()
{
   using blaze::submatrix;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;


   //=====================================================================================
   // Row-major Submatrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major Submatrix multiplication assignment (no aliasing)";

      initialize();

      MT mat1( 64UL, 64UL );
      MT mat2( 64UL, 64UL );
      randomize( mat1, int(randmin), int(randmax) );
      mat2 = mat1;

      ASMT sm1 = submatrix<aligned>  ( mat1_, 16UL, 16UL, 8UL, 8UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 16UL, 16UL, 8UL, 8UL );
      sm1 *= submatrix<aligned>  ( mat1, 16UL, 16UL, 8UL, 8UL );
      sm2 *= submatrix<unaligned>( mat2, 16UL, 16UL, 8UL, 8UL );

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Submatrix multiplication assignment (aliasing)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 16UL, 16UL, 8UL, 8UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 16UL, 16UL, 8UL, 8UL );
      sm1 *= submatrix<aligned>  ( mat1_, 24UL, 16UL, 8UL, 8UL );
      sm2 *= submatrix<unaligned>( mat2_, 24UL, 16UL, 8UL, 8UL );

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix multiplication assignment (mixed type)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 16UL, 16UL, 8UL, 8UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 16UL, 16UL, 8UL, 8UL );

      blaze::DynamicMatrix<short,rowMajor> mat( 8UL, 8UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 16UL, 16UL, 8UL, 8UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 16UL, 16UL, 8UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 8UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major dense matrix multiplication assignment (unaligned/unpadded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 16UL, 16UL, 8UL, 8UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 16UL, 16UL, 8UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[65UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix multiplication assignment (mixed type)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 16UL, 16UL, 8UL, 8UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 16UL, 16UL, 8UL, 8UL );

      blaze::DynamicMatrix<short,columnMajor> mat( 8UL, 8UL );
      randomize( mat, short(randmin), short(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 16UL, 16UL, 8UL, 8UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 16UL, 16UL, 8UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 8UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix multiplication assignment (unaligned/unpadded)";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 16UL, 16UL, 8UL, 8UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 16UL, 16UL, 8UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[65UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix multiplication assignment";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 16UL, 16UL, 8UL, 8UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 16UL, 16UL, 8UL, 8UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 8UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix multiplication assignment";

      initialize();

      ASMT sm1 = submatrix<aligned>  ( mat1_, 16UL, 16UL, 8UL, 8UL );
      USMT sm2 = submatrix<unaligned>( mat2_, 16UL, 16UL, 8UL, 8UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 8UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Submatrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major Submatrix multiplication assignment (no aliasing)";

      initialize();

      OMT mat1( 64UL, 64UL );
      OMT mat2( 64UL, 64UL );
      randomize( mat1, int(randmin), int(randmax) );
      mat2 = mat1;

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 16UL, 8UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 16UL, 8UL, 8UL );
      sm1 *= submatrix<aligned>  ( mat1, 16UL, 16UL, 8UL, 8UL );
      sm2 *= submatrix<unaligned>( mat2, 16UL, 16UL, 8UL, 8UL );

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Submatrix multiplication assignment (aliasing)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 16UL, 8UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 16UL, 8UL, 8UL );
      sm1 *= submatrix<aligned>  ( tmat1_, 16UL, 24UL, 8UL, 8UL );
      sm2 *= submatrix<unaligned>( tmat2_, 16UL, 24UL, 8UL, 8UL );

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix multiplication assignment (mixed type)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 16UL, 8UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 16UL, 8UL, 8UL );

      blaze::DynamicMatrix<short,rowMajor> mat( 8UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 16UL, 8UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 16UL, 8UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 8UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major dense matrix multiplication assignment (unaligned/unpadded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 16UL, 8UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 16UL, 8UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[65UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix multiplication assignment (mixed type)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 16UL, 8UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 16UL, 8UL, 8UL );

      blaze::DynamicMatrix<short,columnMajor> mat( 8UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix multiplication assignment (aligned/padded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 16UL, 8UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 16UL, 8UL, 8UL );

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 128UL ) );
      AlignedPadded mat( memory.get(), 8UL, 8UL, 16UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix multiplication assignment (unaligned/unpadded)";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 16UL, 8UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 16UL, 8UL, 8UL );

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[65UL] );
      UnalignedUnpadded mat( memory.get()+1UL, 8UL, 8UL );
      randomize( mat, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix multiplication assignment";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 16UL, 8UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 16UL, 8UL, 8UL );

      blaze::CompressedMatrix<int,rowMajor> mat( 8UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix multiplication assignment";

      initialize();

      AOSMT sm1 = submatrix<aligned>  ( tmat1_, 16UL, 16UL, 8UL, 8UL );
      UOSMT sm2 = submatrix<unaligned>( tmat2_, 16UL, 16UL, 8UL, 8UL );

      blaze::CompressedMatrix<int,columnMajor> mat( 8UL, 8UL );
      randomize( mat, 30UL, int(randmin), int(randmax) );

      sm1 *= mat;
      sm2 *= mat;

      checkRows   ( sm1, 8UL );
      checkColumns( sm1, 8UL );
      checkRows   ( sm2, 8UL );
      checkColumns( sm2, 8UL );

      if( sm1 != sm2 || mat1_ != mat2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sm1 << "\n"
             << "   Expected result:\n" << sm2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization of all member matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function initializes all member matrices to specific predetermined values.
*/
void DenseAlignedTest::initialize()
{
   // Initializing the row-major dynamic matrices
   randomize( mat1_, int(randmin), int(randmax) );
   mat2_ = mat1_;

   // Initializing the column-major dynamic matrices
   randomize( tmat1_, int(randmin), int(randmax) );
   tmat2_ = tmat1_;
}
//*************************************************************************************************

} // namespace submatrix

} // namespace views

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
   std::cout << "   Running Submatrix dense aligned test (part 1)..." << std::endl;

   try
   {
      RUN_SUBMATRIX_DENSEALIGNED_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during Submatrix dense aligned test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
