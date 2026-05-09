//=================================================================================================
/*!
//  \file src/mathtest/views/subvector/DenseAlignedTest.cpp
//  \brief Source file for the Subvector dense aligned test
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
#include <blaze/math/CompressedVector.h>
#include <blaze/math/CustomVector.h>
#include <blaze/math/Views.h>
#include <blaze/util/Memory.h>
#include <blaze/util/policies/Deallocate.h>
#include <blaze/util/typetraits/AlignmentOf.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>
#include <blazetest/mathtest/views/subvector/DenseAlignedTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace views {

namespace subvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the Subvector dense aligned test.
//
// \exception std::runtime_error Operation error detected.
*/
DenseAlignedTest::DenseAlignedTest()
   : vec1_( 64UL )
   , vec2_( 64UL )
{
   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testDivAssign();
   testCrossAssign();
   testScaling();
   testSubscript();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testIsDefault();
   testIsSame();
   testSubvector();
   testElements();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the Subvector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the Subvector specialization. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testConstructors()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "Subvector constructor";

   initialize();

   const size_t alignment = blaze::AlignmentOf<int>::value;

   for( size_t start=0UL; start<vec1_.size(); start+=alignment ) {
      for( size_t maxsize=0UL; ; maxsize+=alignment )
      {
         const size_t size( blaze::min( maxsize, vec1_.size()-start ) );

         const ASVT sv1 = subvector<aligned>  ( vec1_, start, size );
         const USVT sv2 = subvector<unaligned>( vec2_, start, size );

         if( sv1 != sv2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of dense subvector failed\n"
                << " Details:\n"
                << "   Start = " << start << "\n"
                << "   Size  = " << size  << "\n"
                << "   Subvector:\n" << sv1 << "\n"
                << "   Reference:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( start+maxsize > vec1_.size() ) break;
      }
   }

   try {
      ASVT sv = subvector<aligned>( vec1_, 16UL, 49UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      ASVT sv = subvector<aligned>( vec1_, 80UL, 0UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   if( blaze::AlignmentOf<int>::value > sizeof(int) )
   {
      try {
         ASVT sv = subvector<aligned>( vec1_, 7UL, 16UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of unaligned subvector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the Subvector specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testAssignment()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Homogeneous assignment
   //=====================================================================================

   {
      test_ = "Subvector homogeneous assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      sv1 = 12;
      sv2 = 12;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // List assignment
   //=====================================================================================

   {
      test_ = "Subvector initializer list assignment (complete list)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      sv1 = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 };
      sv2 = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 };

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector initializer list assignment (incomplete list)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      sv1 = { 1, 2, 3 };
      sv2 = { 1, 2, 3 };

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "Subvector copy assignment (no aliasing)";

      initialize();

      VT vec1( 64UL );
      VT vec2( 64UL );
      randomize( vec1, int(randmin), int(randmax) );
      vec2 = vec1;

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      sv1 = subvector<aligned>  ( vec1, 16UL, 21UL );
      sv2 = subvector<unaligned>( vec2, 16UL, 21UL );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector copy assignment (aliasing)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      sv1 = blaze::subvector( vec1_, 32UL, 21UL );
      sv2 = blaze::subvector( vec2_, 32UL, 21UL );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "Subvector dense vector assignment (mixed type)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      blaze::DynamicVector<short,rowVector> vec( 21UL );
      randomize( vec, short(randmin), short(randmax) );

      sv1 = vec;
      sv2 = vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector assignment (aligned/padded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded vec( memory.get(), 21UL, 32UL );
      randomize( vec, int(randmin), int(randmax) );

      sv1 = vec;
      sv2 = vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector assignment (unaligned/unpadded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[22] );
      UnalignedUnpadded vec( memory.get()+1UL, 21UL );
      randomize( vec, int(randmin), int(randmax) );

      sv1 = vec;
      sv2 = vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "Subvector sparse vector assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      blaze::CompressedVector<int,rowVector> vec( 21UL );
      randomize( vec, 6UL, int(randmin), int(randmax) );

      sv1 = vec;
      sv2 = vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testAddAssign()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Subvector addition assignment
   //=====================================================================================

   {
      test_ = "Subvector addition assignment (no aliasing)";

      initialize();

      VT vec1( 64UL );
      VT vec2( 64UL );
      randomize( vec1, int(randmin), int(randmax) );
      vec2 = vec1;

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      sv1 += subvector<aligned>  ( vec1, 16UL, 21UL );
      sv2 += subvector<unaligned>( vec2, 16UL, 21UL );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector addition assignment (aliasing)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      sv1 += subvector<aligned>  ( vec1_, 32UL, 21UL );
      sv2 += subvector<unaligned>( vec2_, 32UL, 21UL );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Subvector dense vector addition assignment (mixed type)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      blaze::DynamicVector<short,rowVector> vec( 21UL );
      randomize( vec, short(randmin), short(randmax) );

      sv1 += vec;
      sv2 += vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector addition assignment (aligned/padded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded vec( memory.get(), 21UL, 32UL );
      randomize( vec, int(randmin), int(randmax) );

      sv1 += vec;
      sv2 += vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector addition assignment (unaligned/unpadded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[22] );
      UnalignedUnpadded vec( memory.get()+1UL, 21UL );
      randomize( vec, int(randmin), int(randmax) );

      sv1 += vec;
      sv2 += vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Subvector sparse vector addition assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      blaze::CompressedVector<int,rowVector> vec( 21UL );
      randomize( vec, 6UL, int(randmin), int(randmax) );

      sv1 += vec;
      sv2 += vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testSubAssign()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Subvector subtraction assignment
   //=====================================================================================

   {
      test_ = "Subvector subtraction assignment (no aliasing)";

      initialize();

      VT vec1( 64UL );
      VT vec2( 64UL );
      randomize( vec1, int(randmin), int(randmax) );
      vec2 = vec1;

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      sv1 -= subvector<aligned>  ( vec1, 32UL, 21UL );
      sv2 -= subvector<unaligned>( vec2, 32UL, 21UL );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector subtraction assignment (aliasing)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      sv1 -= subvector<aligned>  ( vec1_, 32UL, 21UL );
      sv2 -= subvector<unaligned>( vec2_, 32UL, 21UL );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Subvector dense vector subtraction assignment (mixed type)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      blaze::DynamicVector<short,rowVector> vec( 21UL );
      randomize( vec, short(randmin), short(randmax) );

      sv1 -= vec;
      sv2 -= vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector subtraction assignment (aligned/padded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded vec( memory.get(), 21UL, 32UL );
      randomize( vec, int(randmin), int(randmax) );

      sv1 -= vec;
      sv2 -= vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector subtraction assignment (unaligned/unpadded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[22] );
      UnalignedUnpadded vec( memory.get()+1UL, 21UL );
      randomize( vec, int(randmin), int(randmax) );

      sv1 -= vec;
      sv2 -= vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Subvector sparse vector subtraction assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      blaze::CompressedVector<int,rowVector> vec( 21UL );
      randomize( vec, 6UL, int(randmin), int(randmax) );

      sv1 -= vec;
      sv2 -= vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testMultAssign()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Subvector multiplication assignment
   //=====================================================================================

   {
      test_ = "Subvector multiplication assignment (no aliasing)";

      initialize();

      VT vec1( 64UL );
      VT vec2( 64UL );
      randomize( vec1, int(randmin), int(randmax) );
      vec2 = vec1;

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      sv1 *= subvector<aligned>  ( vec1, 32UL, 21UL );
      sv2 *= subvector<unaligned>( vec2, 32UL, 21UL );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector multiplication assignment (aliasing)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      sv1 *= subvector<aligned>  ( vec1_, 32UL, 21UL );
      sv2 *= subvector<unaligned>( vec2_, 32UL, 21UL );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Subvector dense vector multiplication assignment (mixed type)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      blaze::DynamicVector<short,rowVector> vec( 21UL );
      randomize( vec, short(randmin), short(randmax) );

      sv1 *= vec;
      sv2 *= vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector multiplication assignment (aligned/padded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded vec( memory.get(), 21UL, 32UL );
      randomize( vec, int(randmin), int(randmax) );

      sv1 *= vec;
      sv2 *= vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector multiplication assignment (unaligned/unpadded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[22] );
      UnalignedUnpadded vec( memory.get()+1UL, 21UL );
      randomize( vec, int(randmin), int(randmax) );

      sv1 *= vec;
      sv2 *= vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Subvector sparse vector multiplication assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      blaze::CompressedVector<int,rowVector> vec( 21UL );
      randomize( vec, 6UL, int(randmin), int(randmax) );

      sv1 -= vec;
      sv2 -= vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testDivAssign()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Subvector division assignment
   //=====================================================================================

   {
      test_ = "Subvector division assignment (no aliasing)";

      initialize();

      VT vec1( 64UL );
      VT vec2( 64UL );
      randomize( vec1, 1, int(randmax) );
      vec2 = vec1;

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      sv1 /= subvector<aligned>  ( vec1, 32UL, 21UL );
      sv2 /= subvector<unaligned>( vec2, 32UL, 21UL );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector division assignment (aliasing)";

      randomize( vec1_, 1, int(randmax) );
      vec2_ = vec1_;

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      sv1 /= subvector<aligned>  ( vec1_, 32UL, 21UL );
      sv2 /= subvector<unaligned>( vec2_, 32UL, 21UL );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector division assignment
   //=====================================================================================

   {
      test_ = "Subvector dense vector division assignment (mixed type)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      blaze::DynamicVector<short,rowVector> vec( 21UL );
      randomize( vec, short(1), short(randmax) );

      sv1 /= vec;
      sv2 /= vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector division assignment (aligned/padded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded vec( memory.get(), 21UL, 32UL );
      randomize( vec, 1, int(randmax) );

      sv1 /= vec;
      sv2 /= vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector division assignment (unaligned/unpadded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[22] );
      UnalignedUnpadded vec( memory.get()+1UL, 21UL );
      randomize( vec, 1, int(randmax) );

      sv1 /= vec;
      sv2 /= vec;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector cross product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the cross product assignment operators of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testCrossAssign()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowVector;


   //=====================================================================================
   // Subvector cross product assignment
   //=====================================================================================

   {
      test_ = "Subvector cross product assignment (no aliasing)";

      initialize();

      VT vec1( 64UL );
      VT vec2( 64UL );
      randomize( vec1, int(randmin), int(randmax) );
      vec2 = vec1;

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 3UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 3UL );
      sv1 %= subvector<aligned>  ( vec1, 32UL, 3UL );
      sv2 %= subvector<unaligned>( vec2, 32UL, 3UL );

      checkSize( sv1, 3UL );
      checkSize( sv2, 3UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector cross product assignment (aliasing)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 3UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 3UL );
      sv1 %= subvector<aligned>  ( vec1_, 32UL, 3UL );
      sv2 %= subvector<unaligned>( vec2_, 32UL, 3UL );

      checkSize( sv1, 3UL );
      checkSize( sv2, 3UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "Subvector dense vector cross product assignment (mixed type)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 3UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 3UL );

      blaze::DynamicVector<short,rowVector> vec( 3UL );
      randomize( vec, short(randmin), short(randmax) );

      sv1 %= vec;
      sv2 %= vec;

      checkSize( sv1, 3UL );
      checkSize( sv2, 3UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector cross product assignment (aligned/padded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 3UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 3UL );

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec( memory.get(), 3UL, 16UL );
      randomize( vec, int(randmin), int(randmax) );

      sv1 %= vec;
      sv2 %= vec;

      checkSize( sv1, 3UL );
      checkSize( sv2, 3UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Subvector dense vector cross product assignment (unaligned/unpadded)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 3UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 3UL );

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[4] );
      UnalignedUnpadded vec( memory.get()+1UL, 3UL );
      randomize( vec, int(randmin), int(randmax) );

      sv1 %= vec;
      sv2 %= vec;

      checkSize( sv1, 3UL );
      checkSize( sv2, 3UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "Subvector sparse vector cross product assignment";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 3UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 3UL );

      blaze::CompressedVector<int,rowVector> vec( 3UL );
      randomize( vec, 2UL, int(randmin), int(randmax) );

      sv1 %= vec;
      sv2 %= vec;

      checkSize( sv1, 3UL );
      checkSize( sv2, 3UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all Subvector (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testScaling()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   //=====================================================================================
   // Self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "Subvector self-scaling (v*=s)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      sv1 *= 3;
      sv2 *= 3;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v*s)
   //=====================================================================================

   {
      test_ = "Subvector self-scaling (v=v*s)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      sv1 = sv1 * 3;
      sv2 = sv2 * 3;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=s*v)
   //=====================================================================================

   {
      test_ = "Subvector self-scaling (v=s*v)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      sv1 = 3 * sv1;
      sv2 = 3 * sv2;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "Subvector self-scaling (v/=s)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      sv1 /= 0.5;
      sv2 /= 0.5;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v/s)
   //=====================================================================================

   {
      test_ = "Subvector self-scaling (v/=s)";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      sv1 = sv1 / 0.5;
      sv2 = sv2 / 0.5;

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Subvector::scale()
   //=====================================================================================

   {
      test_ = "Subvector::scale()";

      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      // Integral scaling of the subvector in the range [8,23]
      sv1.scale( 3 );
      sv2.scale( 3 );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Integral scale operation of range [8,23] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the subvector in the range [8,23]
      sv1.scale( 0.5 );
      sv2.scale( 0.5 );

      checkSize( sv1, 21UL );
      checkSize( sv2, 21UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Floating point scale operation of range [8,23] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the Subvector specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void DenseAlignedTest::testSubscript()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "Subvector::operator[]";

   initialize();

   ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
   USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

   // Assignment to the element at index 1
   sv1[1] = 9;
   sv2[1] = 9;

   checkSize( sv1, 21UL );
   checkSize( sv2, 21UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 2
   sv1[2] = 0;
   sv2[2] = 0;

   checkSize( sv1, 21UL );
   checkSize( sv2, 21UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 3
   sv1[3] = -8;
   sv2[3] = -8;

   checkSize( sv1, 21UL );
   checkSize( sv2, 21UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Addition assignment to the element at index 0
   sv1[0] += -3;
   sv2[0] += -3;

   checkSize( sv1, 21UL );
   checkSize( sv2, 21UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Subtraction assignment to the element at index 1
   sv1[1] -= 6;
   sv2[1] -= 6;

   checkSize( sv1, 21UL );
   checkSize( sv2, 21UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Multiplication assignment to the element at index 1
   sv1[1] *= 3;
   sv2[1] *= 3;

   checkSize( sv1, 21UL );
   checkSize( sv2, 21UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Division assignment to the element at index 3
   sv1[3] /= 2;
   sv2[3] /= 2;

   checkSize( sv1, 21UL );
   checkSize( sv2, 21UL );

   if( sv1 != sv2 || vec1_ != vec2_ ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Subvector iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the Subvector specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testIterator()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   initialize();

   // Testing the Iterator default constructor
   {
      test_ = "Iterator default constructor";

      ASVT::Iterator it{};

      if( it != ASVT::Iterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing the ConstIterator default constructor
   {
      test_ = "ConstIterator default constructor";

      ASVT::ConstIterator it{};

      if( it != ASVT::ConstIterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing conversion from Iterator to ConstIterator
   {
      test_ = "Iterator/ConstIterator conversion";

      ASVT sv = subvector<aligned>( vec1_, 0UL, 16UL );
      ASVT::ConstIterator it( begin( sv ) );

      if( it == end( sv ) || *it != sv[0] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator conversion detected\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in first half of the vector via Iterator (end-begin)
   {
      test_ = "Iterator subtraction (end-begin)";

      ASVT sv = subvector<aligned>( vec1_, 0UL, 16UL );
      const ptrdiff_t number( end( sv ) - begin( sv ) );

      if( number != 16L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 16\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in first half of the vector via Iterator (begin-end)
   {
      test_ = "Iterator subtraction (begin-end)";

      ASVT sv = subvector<aligned>( vec1_, 0UL, 16UL );
      const ptrdiff_t number( begin( sv ) - end( sv ) );

      if( number != -16L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: -16\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in second half of the vector via ConstIterator (end-begin)
   {
      test_ = "ConstIterator subtraction (end-begin)";

      ASVT sv = subvector<aligned>( vec1_, 16UL, 48UL );
      const ptrdiff_t number( cend( sv ) - cbegin( sv ) );

      if( number != 48L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 48\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in second half of the vector via ConstIterator (begin-end)
   {
      test_ = "ConstIterator subtraction (begin-end)";

      ASVT sv = subvector<aligned>( vec1_, 16UL, 48UL );
      const ptrdiff_t number( cbegin( sv ) - cend( sv ) );

      if( number != -48L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: -48\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      ASVT sv = subvector<aligned>( vec1_, 16UL, 8UL );
      ASVT::ConstIterator it ( cbegin( sv ) );
      ASVT::ConstIterator end( cend( sv ) );

      if( it == end || *it != sv[0] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || *it != sv[1] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      --it;

      if( it == end || *it != sv[0] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it++;

      if( it == end || *it != sv[1] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it--;

      if( it == end || *it != sv[0] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it += 2UL;

      if( it == end || *it != sv[2] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator addition assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it -= 2UL;

      if( it == end || *it != sv[0] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator subtraction assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it + 3UL;

      if( it == end || *it != sv[3] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar addition failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it - 3UL;

      if( it == end || *it != sv[0] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar subtraction failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = 8UL + it;

      if( it != end ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scalar/iterator addition failed\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing assignment via Iterator
   {
      test_ = "Assignment via Iterator";

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      int value = 6;

      ASVT::Iterator it1( begin( sv1 ) );
      USVT::Iterator it2( begin( sv2 ) );

      for( ; it1!=end( sv1 ); ++it1, ++it2 ) {
         *it1 = value;
         *it2 = value;
         ++value;
      }

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing addition assignment via Iterator
   {
      test_ = "Addition assignment via Iterator";

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      int value = 6;

      ASVT::Iterator it1( begin( sv1 ) );
      USVT::Iterator it2( begin( sv2 ) );

      for( ; it1!=end( sv1 ); ++it1, ++it2 ) {
         *it1 += value;
         *it2 += value;
         ++value;
      }

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing subtraction assignment via Iterator
   {
      test_ = "Subtraction assignment via Iterator";

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      int value = 6;

      ASVT::Iterator it1( begin( sv1 ) );
      USVT::Iterator it2( begin( sv2 ) );

      for( ; it1!=end( sv1 ); ++it1, ++it2 ) {
         *it1 -= value;
         *it2 -= value;
         ++value;
      }

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing multiplication assignment via Iterator
   {
      test_ = "Multiplication assignment via Iterator";

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );
      int value = 1;

      ASVT::Iterator it1( begin( sv1 ) );
      USVT::Iterator it2( begin( sv2 ) );

      for( ; it1!=end( sv1 ); ++it1, ++it2 ) {
         *it1 *= value;
         *it2 *= value;
         ++value;
      }

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing division assignment via Iterator
   {
      test_ = "Division assignment via Iterator";

      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
      USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

      ASVT::Iterator it1( begin( sv1 ) );
      USVT::Iterator it2( begin( sv2 ) );

      for( ; it1!=end( sv1 ); ++it1, ++it2 ) {
         *it1 /= 2;
         *it2 /= 2;
      }

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the Subvector specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testNonZeros()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "Subvector::nonZeros()";

   initialize();

   // Initialization check
   ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 21UL );
   USVT sv2 = subvector<unaligned>( vec2_, 16UL, 21UL );

   checkSize( sv1, 21UL );
   checkSize( sv2, 21UL );

   if( sv1.nonZeros() != sv2.nonZeros() ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Changing the number of non-zeros via the dense subvector
   sv1[3] = 0;
   sv2[3] = 0;

   checkSize( sv1, 21UL );
   checkSize( sv2, 21UL );

   if( sv1.nonZeros() != sv2.nonZeros() ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }

   // Changing the number of non-zeros via the dense vector
   vec1_[9UL] = 5;
   vec2_[9UL] = 5;

   checkSize( sv1, 21UL );
   checkSize( sv2, 21UL );

   if( sv1.nonZeros() != sv2.nonZeros() ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << sv1 << "\n"
          << "   Expected result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the Subvector specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testReset()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::reset;


   test_ = "Subvector::reset()";

   // Resetting a single element in the range [0,15]
   {
      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 0UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 0UL, 16UL );
      reset( sv1[4] );
      reset( sv2[4] );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Resetting the range [0,15] (lvalue)
   {
      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 0UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 0UL, 16UL );
      reset( sv1 );
      reset( sv2 );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [0,15] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Resetting the range [16,63] (rvalue)
   {
      initialize();

      reset( subvector<aligned>  ( vec1_, 16UL, 48UL ) );
      reset( subvector<unaligned>( vec2_, 16UL, 48UL ) );

      if( vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [16,63] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1_ << "\n"
             << "   Expected result:\n" << vec2_ << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() function with the Subvector specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() function with the Subvector specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testClear()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::clear;


   test_ = "Subvector::clear()";

   // Clearing a single element in the range [0,15]
   {
      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 0UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 0UL, 16UL );
      clear( sv1[4] );
      clear( sv2[4] );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Clearing the range [0,15] (lvalue)
   {
      initialize();

      ASVT sv1 = subvector<aligned>  ( vec1_, 0UL, 16UL );
      USVT sv2 = subvector<unaligned>( vec2_, 0UL, 16UL );
      clear( sv1 );
      clear( sv2 );

      checkSize( sv1, 16UL );
      checkSize( sv2, 16UL );

      if( sv1 != sv2 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation of range [0,15] failed\n"
             << " Details:\n"
             << "   Result:\n" << sv1 << "\n"
             << "   Expected result:\n" << sv2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Clearing the range [16,63] (rvalue)
   {
      initialize();

      clear( subvector<aligned>  ( vec1_, 16UL, 48UL ) );
      clear( subvector<unaligned>( vec2_, 16UL, 48UL ) );

      if( vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation of range [16,63] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec1_ << "\n"
             << "   Expected result:\n" << vec2_ << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the Subvector specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the Subvector specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testIsDefault()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::isDefault;


   test_ = "isDefault() function";

   initialize();

   // isDefault with default vector
   {
      VT vec( 64UL, 0 );
      ASVT sv = subvector<aligned>( vec, 16UL, 21UL );

      if( isDefault( sv[1] ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Subvector element: " << sv[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( isDefault( sv ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isDefault with non-default vector
   {
      ASVT sv = subvector<aligned>( vec1_, 16UL, 21UL );

      if( isDefault( sv ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Subvector:\n" << sv << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isSame() function with the Subvector specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isSame() function with the Subvector specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testIsSame()
{
   using blaze::subvector;
   using blaze::aligned;


   //=====================================================================================
   // Vector-based tests
   //=====================================================================================

   {
      test_ = "isSame() function (vector-based)";

      // isSame with vector and matching subvector
      {
         ASVT sv = subvector<aligned>( vec1_, 0UL, 64UL );

         if( blaze::isSame( sv, vec1_ ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec1_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( vec1_, sv ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec1_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with vector and non-matching subvector (different size)
      {
         ASVT sv = subvector<aligned>( vec1_, 0UL, 32UL );

         if( blaze::isSame( sv, vec1_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec1_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( vec1_, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec1_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with vector and non-matching subvector (different offset)
      {
         ASVT sv = subvector<aligned>( vec1_, 16UL, 48UL );

         if( blaze::isSame( sv, vec1_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec1_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( vec1_, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec1_ << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching subvectors
      {
         ASVT sv1 = subvector<aligned>( vec1_, 16UL, 32UL );
         ASVT sv2 = subvector<aligned>( vec1_, 16UL, 32UL );

         if( blaze::isSame( sv1, sv2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different size)
      {
         ASVT sv1 = subvector<aligned>( vec1_, 16UL, 32UL );
         ASVT sv2 = subvector<aligned>( vec1_, 16UL, 48UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different offset)
      {
         ASVT sv1 = subvector<aligned>( vec1_, 16UL, 32UL );
         ASVT sv2 = subvector<aligned>( vec1_, 32UL, 32UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-based tests
   //=====================================================================================

   {
      test_ = "isSame() function (row-based)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 64UL, 64UL );
      randomize( mat );

      // isSame with row and matching subvector
      {
         auto r  = blaze::row( mat, 8UL );
         auto sv = subvector<aligned>( r, 0UL, 64UL );

         if( blaze::isSame( sv, r ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( r, sv ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with row and non-matching subvector (different size)
      {
         auto r  = blaze::row( mat, 8UL );
         auto sv = subvector<aligned>( r, 0UL, 32UL );

         if( blaze::isSame( sv, r ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( r, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with row and non-matching subvector (different offset)
      {
         auto r  = blaze::row( mat, 8UL );
         auto sv = subvector<aligned>( r, 16UL, 48UL );

         if( blaze::isSame( sv, r ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( r, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching subvectors
      {
         auto r   = blaze::row( mat, 8UL );
         auto sv1 = subvector<aligned>( r, 0UL, 32UL );
         auto sv2 = subvector<aligned>( r, 0UL, 32UL );

         if( blaze::isSame( sv1, sv2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different size)
      {
         auto r   = blaze::row( mat, 8UL );
         auto sv1 = subvector<aligned>( r, 0UL, 32UL );
         auto sv2 = subvector<aligned>( r, 0UL, 48UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different offset)
      {
         auto r   = blaze::row( mat, 8UL );
         auto sv1 = subvector<aligned>( r,  0UL, 32UL );
         auto sv2 = subvector<aligned>( r, 16UL, 32UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-based tests
   //=====================================================================================

   {
      test_ = "isSame() function (column-based)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 64UL, 64UL );
      randomize( mat );

      // isSame with column and matching subvector
      {
         auto c  = blaze::column( mat, 8UL );
         auto sv = subvector<aligned>( c, 0UL, 64UL );

         if( blaze::isSame( sv, c ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( c, sv ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and non-matching subvector (different size)
      {
         auto c  = blaze::column( mat, 8UL );
         auto sv = subvector<aligned>( c, 0UL, 32UL );

         if( blaze::isSame( sv, c ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( c, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and non-matching subvector (different offset)
      {
         auto c  = blaze::column( mat, 8UL );
         auto sv = subvector<aligned>( c, 16UL, 48UL );

         if( blaze::isSame( sv, c ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( c, sv ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Subvector:\n" << sv << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching subvectors
      {
         auto c   = blaze::column( mat, 8UL );
         auto sv1 = subvector<aligned>( c, 0UL, 32UL );
         auto sv2 = subvector<aligned>( c, 0UL, 32UL );

         if( blaze::isSame( sv1, sv2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different size)
      {
         auto c   = blaze::column( mat, 8UL );
         auto sv1 = subvector<aligned>( c, 0UL, 32UL );
         auto sv2 = subvector<aligned>( c, 0UL, 48UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching subvectors (different offset)
      {
         auto c   = blaze::column( mat, 8UL );
         auto sv1 = subvector<aligned>( c,  0UL, 32UL );
         auto sv2 = subvector<aligned>( c, 16UL, 32UL );

         if( blaze::isSame( sv1, sv2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First subvector:\n" << sv1 << "\n"
                << "   Second subvector:\n" << sv2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c subvector() function with the Subvector specialization.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c subvector() function used with the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testSubvector()
{
   using blaze::subvector;
   using blaze::aligned;
   using blaze::unaligned;


   test_ = "subvector() function";

   initialize();

   {
      ASVT sv1 = subvector<aligned>  ( vec1_, 16UL, 32UL );
      ASVT sv2 = subvector<aligned>  ( sv1  , 16UL, 16UL );
      USVT sv3 = subvector<unaligned>( vec2_, 16UL, 32UL );
      USVT sv4 = subvector<unaligned>( sv3  , 16UL, 16UL );

      if( sv2 != sv4 || vec1_ != vec2_ ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subvector function failed\n"
             << " Details:\n"
             << "   Result:\n" << sv2 << "\n"
             << "   Expected result:\n" << sv4 << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( sv2[1] != sv4[1] ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator access failed\n"
             << " Details:\n"
             << "   Result: " << sv2[1] << "\n"
             << "   Expected result: " << sv4[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( *sv2.begin() != *sv4.begin() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << *sv2.begin() << "\n"
             << "   Expected result: " << *sv4.begin() << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      ASVT sv1 = subvector<aligned>( vec1_, 16UL, 32UL );
      ASVT sv2 = subvector<aligned>( sv1  , 32UL,  8UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      ASVT sv1 = subvector<aligned>( vec1_, 16UL, 32UL );
      ASVT sv2 = subvector<aligned>( sv1  , 16UL, 32UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << sv2 << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c elements() function with the Subvector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c elements() function used with the Subvector
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void DenseAlignedTest::testElements()
{
   using blaze::subvector;
   using blaze::elements;
   using blaze::aligned;
   using blaze::unaligned;


   //=====================================================================================
   // Setup via index_sequence
   //=====================================================================================

   {
      test_ = "elements() function (index_sequence)";

      initialize();

      {
         ASVT sv1 = subvector<aligned>( vec1_, 16UL, 32UL );
         auto e1 = elements( sv1, { 8UL, 16UL } );

         USVT sv2 = subvector<unaligned>( vec2_, 16UL, 32UL );
         auto e2 = elements( sv2, { 8UL, 16UL } );

         if( e1 != e2 || vec1_ != vec2_ ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Elements function failed\n"
                << " Details:\n"
                << "   Result:\n" << e1 << "\n"
                << "   Expected result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( e1[1] != e2[1] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e1[1] << "\n"
                << "   Expected result: " << e2[1] << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e1.begin() != *e2.begin() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e1.begin() << "\n"
                << "   Expected result: " << *e2.begin() << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         ASVT sv = subvector<aligned>( vec1_, 16UL, 32UL );
         auto e = elements( sv, { 8UL, 32UL } );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Setup via std::array
   //=====================================================================================

   {
      test_ = "elements() function (std::array)";

      initialize();

      {
         std::array<int,2UL> indices{ 8UL, 16UL };

         ASVT sv1 = subvector<aligned>( vec1_, 16UL, 32UL );
         auto e1 = elements( sv1, indices );

         USVT sv2 = subvector<unaligned>( vec2_, 16UL, 32UL );
         auto e2 = elements( sv2, indices );

         if( e1 != e2 || vec1_ != vec2_ ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Elements function failed\n"
                << " Details:\n"
                << "   Result:\n" << e1 << "\n"
                << "   Expected result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( e1[1] != e2[1] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e1[1] << "\n"
                << "   Expected result: " << e2[1] << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e1.begin() != *e2.begin() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e1.begin() << "\n"
                << "   Expected result: " << *e2.begin() << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         std::array<int,2UL> indices{ 8UL, 32UL };

         ASVT sv = subvector<aligned>( vec1_, 16UL, 32UL );
         auto e = elements( sv, indices );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Setup via lambda expression
   //=====================================================================================

   {
      test_ = "elements() function (lambda expression)";

      initialize();

      {
         ASVT sv1 = subvector<aligned>( vec1_, 16UL, 32UL );
         auto e1 = elements( sv1, []( size_t i ){ return i*8UL+8UL; }, 2UL );

         USVT sv2 = subvector<unaligned>( vec2_, 16UL, 32UL );
         auto e2 = elements( sv2, []( size_t i ){ return i*8UL+8UL; }, 2UL );

         if( e1 != e2 || vec1_ != vec2_ ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Elements function failed\n"
                << " Details:\n"
                << "   Result:\n" << e1 << "\n"
                << "   Expected result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( e1[1] != e2[1] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e1[1] << "\n"
                << "   Expected result: " << e2[1] << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( *e1.begin() != *e2.begin() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << *e1.begin() << "\n"
                << "   Expected result: " << *e2.begin() << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         ASVT sv = subvector<aligned>( vec1_, 16UL, 32UL );
         auto e = elements( sv, []( size_t i ){ return i*24UL+8UL; }, 2UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization of all member vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function initializes all member vectors to specific predetermined values.
*/
void DenseAlignedTest::initialize()
{
   // Initializing the dynamic row vectors
   randomize( vec1_, int(randmin), int(randmax) );
   vec2_ = vec1_;
}
//*************************************************************************************************

} // namespace subvector

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
   std::cout << "   Running Subvector dense aligned test..." << std::endl;

   try
   {
      RUN_SUBVECTOR_DENSEALIGNED_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during Subvector dense aligned test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
