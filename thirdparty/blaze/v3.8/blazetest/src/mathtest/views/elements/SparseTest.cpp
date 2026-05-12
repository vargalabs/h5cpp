//=================================================================================================
/*!
//  \file src/mathtest/views/elements/SparseTest.cpp
//  \brief Source file for the Elements sparse test
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
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Views.h>
#include <blazetest/mathtest/views/elements/SparseTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace views {

namespace elements {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the Elements sparse test.
//
// \exception std::runtime_error Operation error detected.
*/
SparseTest::SparseTest()
   : vec_( 8UL )
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
   testReserve();
   testSet();
   testInsert();
   testAppend();
   testErase();
   testFind();
   testLowerBound();
   testUpperBound();
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
/*!\brief Test of the Elements constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the Elements specialization. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testConstructors()
{
   using blaze::index_sequence;
   using blaze::initializer_list;


   //=====================================================================================
   // Setup via index_sequence
   //=====================================================================================

   {
      test_ = "Elements constructor (index_sequence)";

      initialize();

      // Setup of a regular element selection
      {
         auto e = blaze::elements( vec_, index_sequence<2,6,4>() );

         if( e.size() != 3UL || e[0] != vec_[2] || e[1] != vec_[6] || e[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds element selection
      try {
         auto e = blaze::elements( vec_, index_sequence<8>() );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of an element selection on a compile-time element selection
      {
         auto e1 = blaze::elements( vec_, index_sequence<2,6,4,3,5>() );
         auto e2 = blaze::elements( e1, index_sequence<1,3,2>() );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a element selection on an explicit element selection
      {
         auto e1 = blaze::elements( vec_, { 2, 6, 4, 3, 5 } );
         auto e2 = blaze::elements( e1, index_sequence<1,3,2>() );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a element selection on an implicit element selection
      {
         const std::array<size_t,5UL> indices{ 2, 6, 4, 3, 5 };
         auto e1 = blaze::elements( vec_, [indices]( size_t i ){ return indices[i]; }, 5UL );
         auto e2 = blaze::elements( e1, index_sequence<1,3,2>() );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Setup via initializer_list
   //=====================================================================================

   {
      test_ = "Elements constructor (initializer_list)";

      initialize();

      // Setup of empty element selection
      {
         std::initializer_list<size_t> indices{};
         auto e = blaze::elements( vec_, indices );

         if( e.size() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of empty element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a regular element selection
      {
         auto e = blaze::elements( vec_, { 2, 6, 4 } );

         if( e.size() != 3UL || e[0] != vec_[2] || e[1] != vec_[6] || e[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds element selection
      try {
         auto e = blaze::elements( vec_, { 8 } );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of an element selection on a compile-time element selection
      {
         auto e1 = blaze::elements( vec_, index_sequence<2,6,4,3,5>() );
         auto e2 = blaze::elements( e1, { 1, 3, 2 } );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a element selection on an explicit element selection
      {
         auto e1 = blaze::elements( vec_, { 2, 6, 4, 3, 5 } );
         auto e2 = blaze::elements( e1, { 1, 3, 2 } );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a element selection on an implicit element selection
      {
         const std::array<size_t,5UL> indices{ 2, 6, 4, 3, 5 };
         auto e1 = blaze::elements( vec_, [indices]( size_t i ){ return indices[i]; }, 5UL );
         auto e2 = blaze::elements( e1, { 1, 3, 2 } );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Setup via std::vector
   //=====================================================================================

   {
      test_ = "Elements constructor (std::vector)";

      initialize();

      // Setup of empty element selection
      {
         const std::vector<size_t> indices;
         auto e = blaze::elements( vec_, indices );

         if( e.size() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of empty element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a regular element selection
      {
         const std::vector<size_t> indices{ 2, 6, 4 };
         auto e = blaze::elements( vec_, indices );

         if( e.size() != 3UL || e[0] != vec_[2] || e[1] != vec_[6] || e[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds element selection
      try {
         const std::vector<size_t> indices{ 8 };
         auto e = blaze::elements( vec_, indices );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of an element selection on a compile-time element selection
      {
         auto e1 = blaze::elements( vec_, index_sequence<2,6,4,3,5>() );

         const std::vector<size_t> indices{ 1, 3, 2 };
         auto e2 = blaze::elements( e1, indices );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a element selection on an explicit element selection
      {
         auto e1 = blaze::elements( vec_, { 2, 6, 4, 3, 5 } );

         const std::vector<size_t> indices{ 1, 3, 2 };
         auto e2 = blaze::elements( e1, indices );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a element selection on an implicit element selection
      {
         const std::array<size_t,5UL> indices1{ 2, 6, 4, 3, 5 };
         auto e1 = blaze::elements( vec_, [indices1]( size_t i ){ return indices1[i]; }, 5UL );

         const std::vector<size_t> indices2{ 1, 3, 2 };
         auto e2 = blaze::elements( e1, indices2 );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Setup via std::array
   //=====================================================================================

   {
      test_ = "Elements constructor (std::array)";

      initialize();

      // Setup of a regular element selection
      {
         const std::array<size_t,3UL> indices{ 2, 6, 4 };
         auto e = blaze::elements( vec_, indices );

         if( e.size() != 3UL || e[0] != vec_[2] || e[1] != vec_[6] || e[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds element selection
      try {
         const std::array<size_t,1UL> indices{ 8 };
         auto e = blaze::elements( vec_, indices );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of an element selection on a compile-time element selection
      {
         auto e1 = blaze::elements( vec_, index_sequence<2,6,4,3,5>() );

         const std::array<size_t,3UL> indices{ 1, 3, 2 };
         auto e2 = blaze::elements( e1, indices );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a element selection on an explicit element selection
      {
         auto e1 = blaze::elements( vec_, { 2, 6, 4, 3, 5 } );

         const std::array<size_t,3UL> indices{ 1, 3, 2 };
         auto e2 = blaze::elements( e1, indices );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a element selection on an implicit element selection
      {
         const std::array<size_t,5UL> indices1{ 2, 6, 4, 3, 5 };
         auto e1 = blaze::elements( vec_, [indices1]( size_t i ){ return indices1[i]; }, 5UL );

         const std::array<size_t,3UL> indices2{ 1, 3, 2 };
         auto e2 = blaze::elements( e1, indices2 );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Setup via lambda expression
   //=====================================================================================

   {
      test_ = "Elements constructor (lambda expression)";

      initialize();

      // Setup of empty element selection
      {
         auto e = blaze::elements( vec_, []( size_t ){ return 0UL; }, 0UL );

         if( e.size() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of empty element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a regular element selection
      {
         const std::array<size_t,3UL> indices{ 2, 6, 4 };
         auto e = blaze::elements( vec_, [indices]( size_t i ){ return indices[i]; }, 3UL );

         if( e.size() != 3UL || e[0] != vec_[2] || e[1] != vec_[6] || e[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds element selection
      try {
         auto e = blaze::elements( vec_, []( size_t ){ return 8UL; }, 1UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds element selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of an element selection on a compile-time element selection
      {
         auto e1 = blaze::elements( vec_, index_sequence<2,6,4,3,5>() );

         const std::array<size_t,3UL> indices{ 1, 3, 2 };
         auto e2 = blaze::elements( e1, [indices]( size_t i ){ return indices[i]; }, 3UL );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a element selection on an explicit element selection
      {
         auto e1 = blaze::elements( vec_, { 2, 6, 4, 3, 5 } );

         const std::array<size_t,3UL> indices{ 1, 3, 2 };
         auto e2 = blaze::elements( e1, [indices]( size_t i ){ return indices[i]; }, 3UL );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a element selection on an implicit element selection
      {
         const std::array<size_t,5UL> indices1{ 2, 6, 4, 3, 5 };
         auto e1 = blaze::elements( vec_, [indices1]( size_t i ){ return indices1[i]; }, 5UL );

         const std::array<size_t,3UL> indices2{ 1, 3, 2 };
         auto e2 = blaze::elements( e1, [indices2]( size_t i ){ return indices2[i]; }, 3UL );

         if( e2.size() != 3UL || e2[0] != vec_[6] || e2[1] != vec_[3] || e2[2] != vec_[4] ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Setup of random in-bounds element selection
   //=====================================================================================

   {
      test_ = "Elements constructor (stress test)";

      initialize();

      for( size_t rep=0UL; rep<100UL; ++rep )
      {
         blaze::DynamicVector<size_t> indices( blaze::rand<size_t>( 1UL, 20UL ) );
         randomize( indices, 0UL, vec_.size()-1UL );
         auto e = blaze::elements( vec_, indices.data(), indices.size() );

         for( size_t i=0UL; i<e.size(); ++i )
         {
            if( e[i] != vec_[indices[i]] ) {
               std::ostringstream oss;
               oss << " Test: " << test_ << "\n"
                   << " Error: Setup of element selection failed\n"
                   << " Details:\n"
                   << "   Indices:\n" << indices << "\n"
                   << "   Element selection:\n" << e << "\n"
                   << "   Vector:\n" << vec_ << "\n";
               throw std::runtime_error( oss.str() );
            }
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Elements assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the Elements specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAssignment()
{
   //=====================================================================================
   // List assignment
   //=====================================================================================

   {
      test_ = "Elements initializer list assignment (complete list)";

      initialize();

      auto e = blaze::elements( vec_, { 0UL, 4UL, 3UL, 7UL } );
      e = { 1, 2, 3, 4 };

      checkSize    ( e   ,  4UL );
      checkNonZeros( e   ,  4UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  6UL );

      if( e[0] != 1 || e[1] != 2 || e[2] != 3 || e[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 1 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 3 ||
          vec_[4] != 2 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 1 1 0 3 2 0 4 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Elements initializer list assignment (incomplete list)";

      initialize();

      auto e = blaze::elements( vec_, { 0UL, 4UL, 3UL, 7UL } );
      e = { 1, 2 };

      checkSize    ( e   ,  4UL );
      checkNonZeros( e   ,  2UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );

      if( e[0] != 1 || e[1] != 2 || e[2] != 0 || e[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 1 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 1 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != 2 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 1 1 0 0 2 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "Elements copy assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      auto e = blaze::elements( vec, { 5UL, 2UL, 7UL } );
      e = blaze::elements( vec_, { 7UL, 3UL, 6UL } );

      checkSize    ( e   ,  3UL );
      checkNonZeros( e   ,  2UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  3UL );

      if( e[0] != 0 || e[1] != -2 || e[2] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 -2 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 0 || vec[1] !=  0 || vec[2] != -2 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != 0 || vec[6] != -8 || vec[7] !=  4 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 -2 0 0 0 -8 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Elements copy assignment (aliasing)";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL } );
      e = blaze::elements( vec_, { 6UL, 5UL, 4UL, 3UL } );

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 4 || e[1] != 0 || e[2] != -3 || e[3] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 4  0 -2 -3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 4 || vec_[2] != 0 || vec_[3] != -3 ||
          vec_[4] != -2 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0  4  0 -3 -2  0  4  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "Dense vector assignment";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 5UL, 6UL } );

      const blaze::DynamicVector<int,blaze::rowVector> vec{ 0, 8, 0, 9 };

      e = vec;

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e != vec ||
          e[0] != 0 || e[1] != 8 || e[2] != 0 || e[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 8 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != 0 || vec_[3] != 8 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 9 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 0 8 -3 0 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "Sparse vector assignment";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 5UL, 6UL } );

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL, 1UL );
      vec[3] = 9;

      e = vec;

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( e != vec ||
          e[0] != 0 || e[1] != 0 || e[2] != 0 || e[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 0 0 9 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 9 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 0 0 -3 0 9 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Elements addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAddAssign()
{
   //=====================================================================================
   // Elements addition assignment
   //=====================================================================================

   {
      test_ = "Elements addition assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      auto e = blaze::elements( vec, { 5UL, 2UL, 7UL } );
      e += blaze::elements( vec_, { 7UL, 3UL, 6UL } );

      checkSize    ( e   ,  3UL );
      checkNonZeros( e   ,  3UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  4UL );

      if( e[0] != 6 || e[1] != -2 || e[2] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 6 -2 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 0 || vec[1] !=  0 || vec[2] != -2 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != 6 || vec[6] != -8 || vec[7] !=  4 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 -2 0 0 6 -8 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Elements addition assignment (aliasing)";

      initialize();

      auto e = blaze::elements( vec_, { 3UL, 4UL, 5UL, 6UL } );
      e += blaze::elements( vec_, { 4UL, 3UL, 2UL, 1UL } );

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != -5 || e[1] != -5 || e[2] != 0 || e[3] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( -5 -5  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != -5 ||
          vec_[4] != -5 || vec_[5] != 0 || vec_[6] != 5 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0  1  0 -5 -5  0  5  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   {
      test_ = "Dense vector addition assignment";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 5UL, 6UL } );

      const blaze::DynamicVector<int,blaze::rowVector> vec{ 0, 8, 0, 9 };

      e += vec;

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 1 || e[1] != 6 || e[2] != 0 || e[3] != 13 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 1 6 0 13 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] !=  0 || vec_[3] != 6 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 13 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 6 -3 0 13 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "Sparse vector addition assignment";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 5UL, 6UL } );

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL, 1UL );
      vec[3] = 9;

      e += vec;

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 1 || e[1] != -2 || e[2] != 0 || e[3] != 13 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 1 -2 0 13 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] !=  0 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 13 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -2 -3 0 13 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Elements subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSubAssign()
{
   //=====================================================================================
   // Elements subtraction assignment
   //=====================================================================================

   {
      test_ = "Elements subtraction assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      auto e = blaze::elements( vec, { 5UL, 2UL, 7UL } );
      e -= blaze::elements( vec_, { 7UL, 3UL, 6UL } );

      checkSize    ( e   ,  3UL );
      checkNonZeros( e   ,  3UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  4UL );

      if( e[0] != 6 || e[1] != 2 || e[2] != -4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 6 2 -4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] != 0 || vec[1] !=  0 || vec[2] !=  2 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != 6 || vec[6] != -8 || vec[7] != -4 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 2 0 0 6 -8 -4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Elements subtraction assignment (aliasing)";

      initialize();

      auto e = blaze::elements( vec_, { 3UL, 4UL, 5UL, 6UL } );
      e -= blaze::elements( vec_, { 4UL, 3UL, 2UL, 1UL } );

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 1 || e[1] != -1 || e[2] != 0 || e[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 1 -1  0  3 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 1 ||
          vec_[4] != -1 || vec_[5] != 0 || vec_[6] != 3 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0  1  0  1 -1  0  3  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Dense vector subtraction assignment";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 5UL, 6UL } );

      const blaze::DynamicVector<int,blaze::rowVector> vec{ 0, 8, 0, 9 };

      e -= vec;

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 1 || e[1] != -10 || e[2] != 0 || e[3] != -5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 1 -10 0 -5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] !=  0 || vec_[3] != -10 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != -5 || vec_[7] !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -10 -3 0 -5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   {
      test_ = "Sparse vector subtraction assignment";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 5UL, 6UL } );

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL, 1UL );
      vec[3] = 9;

      e -= vec;

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 1 || e[1] != -2 || e[2] != 0 || e[3] != -5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 1 -10 0 -5 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] !=  0 || vec_[3] != -2 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != -5 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -10 -3 0 -5 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Elements multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testMultAssign()
{
   //=====================================================================================
   // Elements multiplication assignment
   //=====================================================================================

   {
      test_ = "Elements multiplication assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[5] =  6;
      vec[6] = -8;

      auto e = blaze::elements( vec, { 6UL, 2UL, 5UL } );
      e *= blaze::elements( vec_, { 7UL, 3UL, 6UL } );

      checkSize    ( e   ,  3UL );
      checkNonZeros( e   ,  1UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  1UL );

      if( e[0] != 0 || e[1] != 0 || e[2] != 24 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 0 24 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] !=  0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != 24 || vec[6] != 0 || vec[7] != 0 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 0 24 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Elements multiplication assignment (aliasing)";

      initialize();

      auto e = blaze::elements( vec_, { 3UL, 4UL, 5UL, 6UL } );
      e *= blaze::elements( vec_, { 4UL, 3UL, 2UL, 1UL } );

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 6 || e[1] != 6 || e[2] != 0 || e[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 6  6  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 6 ||
          vec_[4] != 6 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0  1  0  6  6  0  4  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Dense vector multiplication assignment";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 5UL, 6UL } );

      const blaze::DynamicVector<int,blaze::rowVector> vec{ 2, 0, -8, 1 };

      e *= vec;

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e[0] != 2 || e[1] != 0 || e[2] != 0 || e[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 2 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "Sparse vector multiplication assignment";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 5UL, 6UL } );

      blaze::CompressedVector<int,blaze::rowVector> vec( 4UL, 2UL );
      vec[0] = 2;
      vec[3] = 1;

      e *= vec;

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e[0] != 2 || e[1] != 0 || e[2] != 0 || e[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 2 0 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Elements division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testDivAssign()
{
   //=====================================================================================
   // Dense vector division assignment
   //=====================================================================================

   {
      test_ = "Dense vector division assignment";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 5UL, 6UL } );

      const blaze::DynamicVector<int,blaze::rowVector> vec{ 2, -2, 1, -2 };

      e /= vec;

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e[0] != 0 || e[1] != 1 || e[2] != 0 || e[3] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 1 0 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] !=  0 || vec_[3] != 1 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != -2 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 0 1 -3 0 -2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Elements cross product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the cross product assignment operators of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testCrossAssign()
{
   //=====================================================================================
   // Elements cross product assignment
   //=====================================================================================

   {
      test_ = "Elements cross product assignment (no aliasing)";

      initialize();

      VT vec( 10UL, 0 );
      vec[4] =  1;
      vec[6] = -2;
      vec[7] =  4;

      auto e = blaze::elements( vec, { 6UL, 5UL, 4UL } );
      e %= blaze::elements( vec_, { 1UL, 5UL, 3UL } );

      checkSize    ( e   ,  3UL );
      checkNonZeros( e   ,  1UL );
      checkSize    ( vec_,  8UL );
      checkNonZeros( vec_,  4UL );
      checkSize    ( vec , 10UL );
      checkNonZeros( vec ,  2UL );

      if( e[0] != 0 || e[1] != -3 || e[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec[0] !=  0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ||
          vec[5] != -3 || vec[6] != 0 || vec[7] != 4 || vec[8] != 0 || vec[9] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 0 -3 0 4 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Elements cross product assignment (aliasing)";

      initialize();

      auto e = blaze::elements( vec_, { 3UL, 2UL, 1UL } );
      e %= blaze::elements( vec_, { 1UL, 5UL, 3UL } );

      checkSize    ( e   , 3UL );
      checkNonZeros( e   , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e[0] != 0 || e[1] != -3 || e[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != -3 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] !=  4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0  0 -3  0 -3  0  4  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "Dense vector cross product assignment";

      initialize();

      auto e = blaze::elements( vec_, { 3UL, 2UL, 1UL } );

      const blaze::DynamicVector<int,blaze::rowVector> vec{ 1, 0, -2 };

      e %= vec;

      checkSize    ( e   , 3UL );
      checkNonZeros( e   , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e[0] != 0 || e[1] != -3 || e[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != -3 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] !=  4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 -3 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "Sparse vector cross product assignment";

      initialize();

      auto e = blaze::elements( vec_, { 3UL, 2UL, 1UL } );

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] =  1;
      vec[2] = -2;

      e %= vec;

      checkSize    ( e   , 3UL );
      checkNonZeros( e   , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e[0] != 0 || e[1] != -3 || e[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 -3 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != -3 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] !=  4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 -3 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all Elements (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testScaling()
{
   //=====================================================================================
   // Self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "Elements self-scaling (v*=s)";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 2UL } );

      e *= 3;

      checkSize    ( e   , 3UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 3 || e[1] != -6 || e[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 3 -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 0 || vec_[3] != -6 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 0 -6 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v*s)
   //=====================================================================================

   {
      test_ = "Elements self-scaling (v=v*s)";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 2UL } );

      e = e * 3;

      checkSize    ( e   , 3UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 3 || e[1] != -6 || e[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 3 -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 0 || vec_[3] != -6 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 0 -6 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=s*v)
   //=====================================================================================

   {
      test_ = "Elements self-scaling (v=s*v)";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 2UL } );

      e = 3 * e;

      checkSize    ( e   , 3UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 3 || e[1] != -6 || e[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 3 -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 0 || vec_[3] != -6 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 0 -6 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "Elements self-scaling (v/=s)";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 2UL } );

      e /= 0.5;

      checkSize    ( e   , 3UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 2 || e[1] != -4 || e[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 3 -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != -4 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 -4 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v/s)
   //=====================================================================================

   {
      test_ = "Elements self-scaling (v=v/s)";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 2UL } );

      e = e / 0.5;

      checkSize    ( e   , 3UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 2 || e[1] != -4 || e[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 3 -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 2 || vec_[2] != 0 || vec_[3] != -4 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 2 0 -4 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Elements::scale()
   //=====================================================================================

   {
      test_ = "Elements::scale()";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 3UL, 2UL } );

      // Integral scaling of the element selection in the range [1,4]
      e.scale( 3 );

      checkSize    ( e   , 3UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 3 || e[1] != -6 || e[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 3 -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 3 || vec_[2] != 0 || vec_[3] != -6 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 3 0 -6 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the element selection in the range [1,4]
      e.scale( 0.5 );

      checkSize    ( e   , 3UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 4UL );

      if( e[0] != 1 || e[1] != -3 || e[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 1 -3  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != -3 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -3 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Elements subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the Elements specialization. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void SparseTest::testSubscript()
{
   test_ = "Elements::operator[]";

   initialize();

   auto e = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL } );

   // Assignment to the element at index 1
   e[1] = 9;

   checkSize    ( e   , 4UL );
   checkNonZeros( e   , 4UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 5UL );

   if( e[0] != 1 || e[1] != 9 || e[2] != -2 || e[3] != -3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( 1 9 -2 -3 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 9 || vec_[3] != -2 ||
       vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] !=  0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 1 9 -2 -3 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 2
   e[2] = 0;

   checkSize    ( e   , 4UL );
   checkNonZeros( e   , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( e[0] != 1 || e[1] != 9 || e[2] != 0 || e[3] != -3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( 1 9 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 9 || vec_[3] != 0 ||
       vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 1 9 0 -3 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 8
   e[3] = -8;

   checkSize    ( e   , 4UL );
   checkNonZeros( e   , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( e[0] != 1 || e[1] != 9 || e[2] != 0 || e[3] != -8 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( 1 9 0 -8 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 9 || vec_[3] != 0 ||
       vec_[4] != -8 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 1 9 0 -8 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Addition assignment to the element at index 0
   e[0] += -3;

   checkSize    ( e   , 4UL );
   checkNonZeros( e   , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( e[0] != -2 || e[1] != 9 || e[2] != 0 || e[3] != -8 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( -2 9 0 -8 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != -2 || vec_[2] != 9 || vec_[3] != 0 ||
       vec_[4] != -8 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 -2 9 0 -8 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Subtraction assignment to the element at index 1
   e[1] -= 6;

   checkSize    ( e   , 4UL );
   checkNonZeros( e   , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( e[0] != -2 || e[1] != 3 || e[2] != 0 || e[3] != -8 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( -2 3 0 -8 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != -2 || vec_[2] != 3 || vec_[3] != 0 ||
       vec_[4] != -8 || vec_[5] !=  0 || vec_[6] != 4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 -2 3 0 -8 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Multiplication assignment to the element at index 1
   e[1] *= -3;

   checkSize    ( e   , 4UL );
   checkNonZeros( e   , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( e[0] != -2 || e[1] != -9 || e[2] != 0 || e[3] != -8 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( -2 -9 0 -8 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != -2 || vec_[2] != -9 || vec_[3] != 0 ||
       vec_[4] != -8 || vec_[5] !=  0 || vec_[6] !=  4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 -2 -9 0 -8 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Multiplication assignment to the element at index 3
   e[3] /= 2;

   checkSize    ( e   , 4UL );
   checkNonZeros( e   , 3UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( e[0] != -2 || e[1] != -9 || e[2] != 0 || e[3] != -4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( -2 -9 0 -4 )\n";
      throw std::runtime_error( oss.str() );
   }

   if( vec_[0] !=  0 || vec_[1] != -2 || vec_[2] != -9 || vec_[3] != 0 ||
       vec_[4] != -4 || vec_[5] !=  0 || vec_[6] !=  4 || vec_[7] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec_ << "\n"
          << "   Expected result:\n( 0 -2 -9 0 -4 0 4 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Elements iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the Elements specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testIterator()
{
   initialize();

   // Testing the Iterator default constructor
   {
      test_ = "Iterator default constructor";

      ET::Iterator it{};

      if( it != ET::Iterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing the ConstIterator default constructor
   {
      test_ = "ConstIterator default constructor";

      ET::ConstIterator it{};

      if( it != ET::ConstIterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing conversion from Iterator to ConstIterator
   {
      test_ = "Iterator/ConstIterator conversion";

      auto e = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL } );
      auto it( begin( e ) );

      if( it == end( e ) || it->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator conversion detected\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in first half of the vector via Iterator (end-begin)
   {
      test_ = "Iterator subtraction (end-begin)";

      auto e = blaze::elements( vec_, { 0UL, 1UL, 2UL, 3UL } );
      const ptrdiff_t number( end( e ) - begin( e ) );

      if( number != 2L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 2\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in second half of the vector via ConstIterator (end-begin)
   {
      test_ = "ConstIterator subtraction (end-begin)";

      auto e = blaze::elements( vec_, { 4UL, 5UL, 6UL, 7UL } );
      const ptrdiff_t number( cend( e ) - cbegin( e ) );

      if( number != 2L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 2\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      auto e = blaze::elements( vec_, { 1UL, 2UL, 3UL } );
      auto it ( cbegin( e ) );
      auto end( cend( e ) );

      if( it == end || it->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || it->value() != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it++;

      if( it != end ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-increment failed\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing assignment via Iterator
   {
      test_ = "Assignment via Iterator";

      auto e = blaze::elements( vec_, { 2UL, 3UL, 4UL, 5UL } );
      int value = 6;

      for( auto it=begin( e ); it!=end( e ); ++it ) {
         *it = value++;
      }

      if( e[0] != 0 || e[1] != 6 || e[2] != 7 || e[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 6 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 6 ||
          vec_[4] != 7 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 6 7 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing addition assignment via Iterator
   {
      test_ = "Addition assignment via Iterator";

      auto e = blaze::elements( vec_, { 2UL, 3UL, 4UL, 5UL } );
      int value = 2;

      for( auto it=begin( e ); it!=end( e ); ++it ) {
         *it += value++;
      }

      if( e[0] != 0 || e[1] != 8 || e[2] != 10 || e[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 8 10 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 8 ||
          vec_[4] != 10 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 8 10 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing subtraction assignment via Iterator
   {
      test_ = "Subtraction assignment via Iterator";

      auto e = blaze::elements( vec_, { 2UL, 3UL, 4UL, 5UL } );
      int value = 2;

      for( auto it=begin( e ); it!=end( e ); ++it ) {
         *it -= value++;
      }

      if( e[0] != 0 || e[1] != 6 || e[2] != 7 || e[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 6 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 6 ||
          vec_[4] != 7 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 6 7 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing multiplication assignment via Iterator
   {
      test_ = "Multiplication assignment via Iterator";

      auto e = blaze::elements( vec_, { 2UL, 3UL, 4UL, 5UL } );
      int value = 1;

      for( auto it=begin( e ); it!=end( e ); ++it ) {
         *it *= value++;
      }

      if( e[0] != 0 || e[1] != 6 || e[2] != 14 || e[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 6 14 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 6 ||
          vec_[4] != 14 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 6 14 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing division assignment via Iterator
   {
      test_ = "Division assignment via Iterator";

      auto e = blaze::elements( vec_, { 2UL, 3UL, 4UL, 5UL } );

      for( auto it=begin( e ); it!=end( e ); ++it ) {
         *it /= 2;
      }

      if( e[0] != 0 || e[1] != 3 || e[2] != 7 || e[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 3 7 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != 3 ||
          vec_[4] != 7 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 3 7 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testNonZeros()
{
   test_ = "Elements::nonZeros()";

   initialize();

   // Initialization check
   auto e = blaze::elements( vec_, { 3UL, 2UL, 1UL, 0UL } );

   checkSize    ( e   , 4UL );
   checkNonZeros( e   , 2UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( e[0] != -2 || e[1] != 0 || e[2] != 1 || e[3] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Initialization failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( -2 0 1 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Changing the number of non-zeros via the dense element selection
   e[0] = 0;

   checkSize    ( e   , 4UL );
   checkNonZeros( e   , 1UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 3UL );

   if( e[0] != 0 || e[1] != 0 || e[2] != 1 || e[3] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( 0 0 1 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Changing the number of non-zeros via the dense vector
   vec_[2UL] = 5;

   checkSize    ( e   , 4UL );
   checkNonZeros( e   , 2UL );
   checkSize    ( vec_, 8UL );
   checkNonZeros( vec_, 4UL );

   if( e[0] != 0 || e[1] != 5 || e[2] != 1 || e[3] != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( 0 5 1 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testReset()
{
   test_ = "Elements::reset()";

   using blaze::reset;

   // Resetting a single element of the range [1,6]
   {
      initialize();

      auto e = blaze::elements( vec_, { 6UL, 3UL, 2UL, 5UL, 4UL, 1UL } );
      reset( e[1] );

      checkSize    ( e   , 6UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e[0] != 4 || e[1] != 0 || e[2] != 0 || e[3] != 0 || e[4] != -3 || e[5] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 4 0 0 0 -3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Resetting the range [0,3] (lvalue)
   {
      initialize();

      auto e = blaze::elements( vec_, { 3UL, 2UL, 1UL, 0UL } );
      reset( e );

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 0UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( e[0] != 0 || e[1] != 0 || e[2] != 0 || e[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [0,3] failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [0,3] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Resetting the range [4,7] (rvalue)
   {
      initialize();

      reset( blaze::elements( vec_, { 4UL, 5UL, 6UL, 7UL } ) );

      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != -2 ||
          vec_[4] != 0 || vec_[5] != 0 || vec_[6] != 0 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation of range [4,7] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -2 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() function with the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() function of the Elements specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testClear()
{
   test_ = "clear() function";

   using blaze::clear;

   // Clearing a single element of the range [1,6]
   {
      initialize();

      auto e = blaze::elements( vec_, { 6UL, 3UL, 2UL, 5UL, 4UL, 1UL } );
      clear( e[1] );

      checkSize    ( e   , 6UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e[0] != 4 || e[1] != 0 || e[2] != 0 || e[3] != 0 || e[4] != -3 || e[5] != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 4 0 0 0 -3 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Clearing the range [0,3] (lvalue)
   {
      initialize();

      auto e = blaze::elements( vec_, { 3UL, 2UL, 1UL, 0UL } );
      clear( e );

      checkSize    ( e   , 4UL );
      checkNonZeros( e   , 0UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( e[0] != 0 || e[1] != 0 || e[2] != 0 || e[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation of range [0,3] failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( vec_[0] !=  0 || vec_[1] != 0 || vec_[2] != 0 || vec_[3] != 0 ||
          vec_[4] != -3 || vec_[5] != 0 || vec_[6] != 4 || vec_[7] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation of range [0,3] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 0 0 0 -3 0 4 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Clearing the range [4,7] (rvalue)
   {
      initialize();

      clear( blaze::elements( vec_, { 4UL, 5UL, 6UL, 7UL } ) );

      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( vec_[0] != 0 || vec_[1] != 1 || vec_[2] != 0 || vec_[3] != -2 ||
          vec_[4] != 0 || vec_[5] != 0 || vec_[6] != 0 || vec_[7] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation of range [4,7] failed\n"
             << " Details:\n"
             << "   Result:\n" << vec_ << "\n"
             << "   Expected result:\n( 0 1 0 -2 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reserve() member function of the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reserve() member function of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testReserve()
{
   test_ = "Elements::reserve()";

   VT vec( 10UL );

   auto e = blaze::elements( vec, { 2UL, 5UL, 4UL, 3UL } );

   // Increasing the capacity of the vector
   e.reserve( 10UL );

   checkSize    ( e,  4UL );
   checkCapacity( e, 10UL );
   checkNonZeros( e,  0UL );

   // Further increasing the capacity of the vector
   e.reserve( 20UL );

   checkSize    ( e,  4UL );
   checkCapacity( e, 20UL );
   checkNonZeros( e,  0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c set() member function of the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c set() member function of the Elements specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSet()
{
   test_ = "Elements::set()";

   initialize();

   auto e = blaze::elements( vec_, { 0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 6UL, 7UL } );

   // Setting a non-zero element at the end of the element selection
   {
      auto pos = e.set( 7UL, 9 );

      checkSize    ( e   , 8UL );
      checkNonZeros( e   , 5UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( pos->value() != 9 || pos->index() != 7UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 7\n";
         throw std::runtime_error( oss.str() );
      }

      if( e[0] !=  0 || e[1] != 1 || e[2] != 0 || e[3] != -2 ||
          e[4] != -3 || e[5] != 0 || e[6] != 4 || e[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 2 0 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Setting a non-zero element at the beginning of the element selection
   {
      auto pos = e.set( 0UL, 9 );

      checkSize    ( e   , 8UL );
      checkNonZeros( e   , 6UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 6UL );

      if( pos->value() != 9 || pos->index() != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 0\n";
         throw std::runtime_error( oss.str() );
      }

      if( e[0] !=  9 || e[1] != 1 || e[2] != 0 || e[3] != -2 ||
          e[4] != -3 || e[5] != 0 || e[6] != 4 || e[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 9 2 0 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Setting a non-zero element at the center of the element selection
   {
      auto pos = e.set( 2UL, 9 );

      checkSize    ( e   , 8UL );
      checkNonZeros( e   , 7UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 7UL );

      if( pos->value() != 9 || pos->index() != 2UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 2\n";
         throw std::runtime_error( oss.str() );
      }

      if( e[0] !=  9 || e[1] != 1 || e[2] != 9 || e[3] != -2 ||
          e[4] != -3 || e[5] != 0 || e[6] != 4 || e[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 9 2 9 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Setting an already existing element
   {
      auto pos = e.set( 3UL, 9 );

      checkSize    ( e   , 8UL );
      checkNonZeros( e   , 7UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 7UL );

      if( pos->value() != 9 || pos->index() != 3UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 3\n";
         throw std::runtime_error( oss.str() );
      }

      if( e[0] !=  9 || e[1] != 1 || e[2] != 9 || e[3] != 9 ||
          e[4] != -3 || e[5] != 0 || e[6] != 4 || e[7] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 9 2 9 9 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c insert() member function of the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c insert() member function of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testInsert()
{
   test_ = "Elements::insert()";

   initialize();

   auto e = blaze::elements( vec_, { 0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 6UL, 7UL } );

   // Inserting a non-zero element at the end of the element selection
   {
      auto pos = e.insert( 7UL, 9 );

      checkSize    ( e   , 8UL );
      checkNonZeros( e   , 5UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 5UL );

      if( pos->value() != 9 || pos->index() != 7UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 7\n";
         throw std::runtime_error( oss.str() );
      }

      if( e[0] !=  0 || e[1] != 1 || e[2] != 0 || e[3] != -2 ||
          e[4] != -3 || e[5] != 0 || e[6] != 4 || e[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 2 0 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a non-zero element at the beginning of the element selection
   {
      auto pos = e.insert( 0UL, 9 );

      checkSize    ( e   , 8UL );
      checkNonZeros( e   , 6UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 6UL );

      if( pos->value() != 9 || pos->index() != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 0\n";
         throw std::runtime_error( oss.str() );
      }

      if( e[0] !=  9 || e[1] != 1 || e[2] != 0 || e[3] != -2 ||
          e[4] != -3 || e[5] != 0 || e[6] != 4 || e[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 9 2 0 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Inserting a non-zero element at the center of the element selection
   {
      auto pos = e.insert( 2UL, 9 );

      checkSize    ( e   , 8UL );
      checkNonZeros( e   , 7UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 7UL );

      if( pos->value() != 9 || pos->index() != 2UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid iterator returned\n"
             << " Details:\n"
             << "   Value: " << pos->value() << "\n"
             << "   Index: " << pos->index() << "\n"
             << "   Expected value: 9\n"
             << "   Expected index: 2\n";
         throw std::runtime_error( oss.str() );
      }

      if( e[0] !=  9 || e[1] != 1 || e[2] != 9 || e[3] != -2 ||
          e[4] != -3 || e[5] != 0 || e[6] != 4 || e[7] !=  9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Inserting a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 9 2 9 -2 -3 0 4 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Trying to insert an already existing element
   try {
      e.insert( 3UL, 9 );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Inserting an existing element succeeded\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( 9 2 0 9 -3 0 4 9 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c append() member function of the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c append() member function of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testAppend()
{
   test_ = "Elements::append()";

   VT vec( 10UL );

   auto e = blaze::elements( vec, { 2UL, 3UL, 4UL, 5UL } );
   e.reserve( 4UL );

   // Appending one non-zero element
   e.append( 0UL, 1 );

   checkSize    ( e  , 4UL );
   checkCapacity( e  , 4UL );
   checkNonZeros( e  , 1UL );
   checkNonZeros( vec, 1UL );

   if( e[0] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Append operation failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( 0 0 1 0 0 0 0 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Appending three more non-zero elements
   e.append( 1UL, 2 );
   e.append( 2UL, 3 );
   e.append( 3UL, 4 );

   checkSize    ( e  , 4UL );
   checkCapacity( e  , 4UL );
   checkNonZeros( e  , 4UL );
   checkNonZeros( vec, 4UL );

   if( e[0] != 1 || e[1] != 2 || e[2] != 3 || e[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Append operation failed\n"
          << " Details:\n"
          << "   Result:\n" << e << "\n"
          << "   Expected result:\n( 0 0 1 2 3 4 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c erase() member function of the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c erase() member function of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testErase()
{
   //=====================================================================================
   // Index-based erase function
   //=====================================================================================

   {
      test_ = "Elements::erase( size_t )";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL, 5UL, 6UL } );

      // Erasing the non-zero element at the end of the element selection
      e.erase( 5UL );

      checkSize    ( e   , 6UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e[0] !=  1 || e[1] != 0 || e[2] != -2 ||
          e[3] != -3 || e[4] != 0 || e[5] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 1 0 -2 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the element selection
      e.erase( 0UL );

      checkSize    ( e   , 6UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( e[0] !=  0 || e[1] != 0 || e[2] != -2 ||
          e[3] != -3 || e[4] != 0 || e[5] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 0 -2 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Erasing the non-zero element at the beginning of the element selection
      e.erase( 2UL );

      checkSize    ( e   , 6UL );
      checkNonZeros( e   , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 1UL );

      if( e[0] !=  0 || e[1] != 0 || e[2] != 0 ||
          e[3] != -3 || e[4] != 0 || e[5] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a non-zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase an already erased element
      e.erase( 2UL );

      checkSize    ( e   , 6UL );
      checkNonZeros( e   , 1UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 1UL );

      if( e[0] !=  0 || e[1] != 0 || e[2] != 0 ||
          e[3] != -3 || e[4] != 0 || e[5] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a zero element failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Iterator-based erase function
   //=====================================================================================

   {
      test_ = "Elements::erase( Iterator )";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL, 5UL, 6UL } );

      // Erasing the non-zero element at the end of the element selection
      {
         auto pos = e.erase( e.find( 5UL ) );

         checkSize    ( e   , 6UL );
         checkNonZeros( e   , 3UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 3UL );

         if( pos != e.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( e[0] !=  1 || e[1] != 0 || e[2] != -2 ||
             e[3] != -3 || e[4] != 0 || e[5] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( 1 0 -2 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the element selection
      {
         auto pos = e.erase( e.find( 0UL ) );

         checkSize    ( e   , 6UL );
         checkNonZeros( e   , 2UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 2UL );

         if( pos->value() != -2 || pos->index() != 2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -2\n"
                << "   Expected index:  2\n";
            throw std::runtime_error( oss.str() );
         }

         if( e[0] !=  0 || e[1] != 0 || e[2] != -2 ||
             e[3] != -3 || e[4] != 0 || e[5] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( 0 0 -2 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the non-zero element at the beginning of the element selection
      {
         auto pos = e.erase( e.find( 2UL ) );

         checkSize    ( e   , 6UL );
         checkNonZeros( e   , 1UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 1UL );

         if( pos->value() != -3 || pos->index() != 3 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -3\n"
                << "   Expected index:  3\n";
            throw std::runtime_error( oss.str() );
         }

         if( e[0] !=  0 || e[1] != 0 || e[2] != 0 ||
             e[3] != -3 || e[4] != 0 || e[5] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a non-zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an already erased element
      {
         auto pos = e.erase( e.find( 2UL ) );

         checkSize    ( e   , 6UL );
         checkNonZeros( e   , 1UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 1UL );

         if( pos != e.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( e[0] !=  0 || e[1] != 0 || e[2] != 0 ||
             e[3] != -3 || e[4] != 0 || e[5] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a zero element failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( 0 0 0 -3 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Iterator-range-based erase function
   //=====================================================================================

   {
      test_ = "Elements::erase( Iterator, Iterator )";

      // Erasing the entire vector
      {
         initialize();

         auto e = blaze::elements( vec_, { 0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 6UL, 7UL } );

         auto pos = e.erase( e.begin(), e.end() );

         checkSize    ( e   , 8UL );
         checkNonZeros( e   , 0UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 0UL );

         if( pos != e.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( e[0] != 0 || e[1] != 0 || e[2] != 0 || e[3] != 0 ||
             e[4] != 0 || e[5] != 0 || e[6] != 0 || e[7] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing the complete element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( 0 0 0 0 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the first half of the vector
      {
         initialize();

         auto e = blaze::elements( vec_, { 0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 6UL, 7UL } );

         auto pos = e.erase( e.begin(), e.find( 4UL ) );

         checkSize    ( e   , 8UL );
         checkNonZeros( e   , 2UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 2UL );

         if( pos->value() != -3 || pos->index() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Value: " << pos->value() << "\n"
                << "   Index: " << pos->index() << "\n"
                << "   Expected value: -3\n"
                << "   Expected index:  4\n";
            throw std::runtime_error( oss.str() );
         }

         if( e[0] !=  0 || e[1] != 0 || e[2] != 0 || e[3] != 0 ||
             e[4] != -3 || e[5] != 0 || e[6] != 4 || e[7] != 0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( 0 0 0 0 -3 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Erasing the second half of the vector
      {
         initialize();

         auto e = blaze::elements( vec_, { 0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 6UL, 7UL } );

         auto pos = e.erase( e.find( 4UL ), e.end() );

         checkSize    ( e   , 8UL );
         checkNonZeros( e   , 2UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 2UL );

         if( pos != e.end() ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( e[0] != 0 || e[1] != 1 || e[2] != 0 || e[3] != -2 ||
             e[4] != 0 || e[5] != 0 || e[6] != 0 || e[7] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing a partial element selection failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( 0 1 0 -2 0 0 0 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to erase an empty range
      {
         initialize();

         auto e = blaze::elements( vec_, { 0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 6UL, 7UL } );

         auto pos = e.erase( e.find( 1UL ), e.find( 1UL ) );

         checkSize    ( e   , 8UL );
         checkNonZeros( e   , 4UL );
         checkSize    ( vec_, 8UL );
         checkNonZeros( vec_, 4UL );

         if( pos != e.find( 1UL ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid iterator returned\n"
                << " Details:\n"
                << "   Expected result: the end() iterator\n";
            throw std::runtime_error( oss.str() );
         }

         if( e[0] !=  0 || e[1] != 1 || e[2] != 0 || e[3] != -2 ||
             e[4] != -3 || e[5] != 0 || e[6] != 4 || e[7] !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Erasing an empty range failed\n"
                << " Details:\n"
                << "   Result:\n" << e << "\n"
                << "   Expected result:\n( 0 1 0 -2 -3 0 4 0 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   //  erase() function with predicate
   //=====================================================================================

   {
      test_ = "Elements::erase( Predicate )";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL, 5UL, 6UL } );

      // Erasing a selection of elements
      e.erase( []( int value ){ return value == 1 || value == 4; } );

      checkSize    ( e   , 6UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( e[0] !=  0 || e[1] != 0 || e[2] != -2 ||
          e[3] != -3 || e[4] != 0 || e[5] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 0 -2 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase all elements with value 1
      e.erase( []( int value ){ return value == 1; } );

      checkSize    ( e   , 6UL );
      checkNonZeros( e   , 2UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 2UL );

      if( e[0] !=  0 || e[1] != 0 || e[2] != -2 ||
          e[3] != -3 || e[4] != 0 || e[5] !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing all elements with value 1 failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 0 -2 -3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Iterator-range-based erase() function with predicate
   //=====================================================================================

   {
      test_ = "Elements::erase( Iterator, Iterator, Predicate )";

      initialize();

      auto e = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL, 5UL, 6UL } );

      // Erasing a selection of elements
      e.erase( e.begin(), e.find( 3UL ), []( int value ){ return value == 1; } );

      checkSize    ( e   , 6UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e[0] !=  0 || e[1] != 0 || e[2] != -2 ||
          e[3] != -3 || e[4] != 0 || e[5] !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing a selection of elements failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 0 -2 -3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Trying to erase from an empty range
      e.erase( e.begin(), e.begin(), []( int value ){ return value == 1; } );

      checkSize    ( e   , 6UL );
      checkNonZeros( e   , 3UL );
      checkSize    ( vec_, 8UL );
      checkNonZeros( vec_, 3UL );

      if( e[0] !=  0 || e[1] != 0 || e[2] != -2 ||
          e[3] != -3 || e[4] != 0 || e[5] !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Erasing from an empty range failed\n"
             << " Details:\n"
             << "   Result:\n" << e << "\n"
             << "   Expected result:\n( 0 0 -2 -3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c find() member function of the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c find() member function of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testFind()
{
   test_ = "Elements::find()";

   initialize();

   auto e = blaze::elements( vec_, { 1UL, 5UL, 2UL, 4UL, 3UL } );

   // Searching for the first element
   {
      auto pos = e.find( 0UL );

      if( pos == e.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 0 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for the second element
   {
      auto pos = e.find( 3UL );

      if( pos == e.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Element could not be found\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 3 || pos->value() != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 3\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = -3\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Searching for a non-existing non-zero element
   {
      auto pos = e.find( 1UL );

      if( pos != e.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Non-existing element could be found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c lowerBound() member function of the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c lowerBound() member function of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testLowerBound()
{
   test_ = "Elements::lowerBound()";

   initialize();

   auto e = blaze::elements( vec_, { 2UL, 1UL, 0UL } );

   // Determining the lower bound for index 0
   {
      auto pos = e.lowerBound( 0UL );

      if( pos == e.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 1 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 1
   {
      auto pos = e.lowerBound( 1UL );

      if( pos == e.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 1 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the lower bound for index 2
   {
      auto pos = e.lowerBound( 2UL );

      if( pos != e.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Lower bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c upperBound() member function of the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c upperBound() member function of the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testUpperBound()
{
   test_ = "Elements::upperBound()";

   initialize();

   auto e = blaze::elements( vec_, { 2UL, 1UL, 0UL } );

   // Determining the upper bound for index 0
   {
      auto pos = e.upperBound( 0UL );

      if( pos == e.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 0\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
      else if( pos->index() != 1 || pos->value() != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Wrong element found\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Found index    = " << pos->index() << "\n"
             << "   Expected value = 1\n"
             << "   Value at index = " << pos->value() << "\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 1
   {
      auto pos = e.upperBound( 1UL );

      if( pos != e.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 1\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Determining the upper bound for index 2
   {
      auto pos = e.upperBound( 2UL );

      if( pos != e.end() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Upper bound could not be determined\n"
             << " Details:\n"
             << "   Required index = 2\n"
             << "   Current element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the Elements specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testIsDefault()
{
   test_ = "isDefault() function";

   using blaze::isDefault;

   initialize();

   // isDefault with default vector
   {
      VT vec( 8UL );
      auto e = blaze::elements( vec, { 5UL, 4UL, 6UL, 2UL, 3UL } );

      if( isDefault( e[1] ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Element: " << e[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( isDefault( e ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isDefault with non-default vector
   {
      auto e = blaze::elements( vec_, { 5UL, 4UL, 6UL, 2UL, 3UL } );

      if( isDefault( e[1] ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Element: " << e[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( isDefault( e ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Element selection:\n" << e << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isSame() function with the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isSame() function with the Elements specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testIsSame()
{
   //=====================================================================================
   // Vector-based tests
   //=====================================================================================

   {
      test_ = "isSame() function (vector-based)";

      // isSame with vector and matching element selection
      {
         auto e = blaze::elements( vec_, { 0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 6UL, 7UL } );

         if( blaze::isSame( e, vec_ ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( vec_, e ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with vector and non-matching element selection (different size)
      {
         auto e = blaze::elements( vec_, { 0UL, 1UL, 2UL, 3UL, 4UL, 5UL, 6UL } );

         if( blaze::isSame( e, vec_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( vec_, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with vector and non-matching element selection (different order)
      {
         auto e = blaze::elements( vec_, { 0UL, 1UL, 3UL, 2UL, 4UL, 5UL, 6UL, 7UL } );

         if( blaze::isSame( e, vec_ ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( vec_, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Vector:\n" << vec_ << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with subvector and matching element selection
      {
         auto e = blaze::elements( vec_, { 2UL, 3UL, 4UL } );
         auto s = blaze::subvector( vec_, 2UL, 3UL );

         if( blaze::isSame( e, s ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subvector:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( s, e ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subvector:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with subvector and non-matching element selection (different size)
      {
         auto e = blaze::elements( vec_, { 2UL, 3UL, 4UL } );
         auto s = blaze::subvector( vec_, 2UL, 4UL );

         if( blaze::isSame( e, s ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subvector:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( s, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subvector:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with subvector and non-matching element selection (different order)
      {
         auto e = blaze::elements( vec_, { 2UL, 4UL, 3UL } );
         auto s = blaze::subvector( vec_, 2UL, 3UL );

         if( blaze::isSame( e, s ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subvector:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( s, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subvector:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching element selections
      {
         auto e1 = blaze::elements( vec_, { 5UL, 3UL, 1UL } );
         auto e2 = blaze::elements( vec_, { 5UL, 3UL, 1UL } );

         if( blaze::isSame( e1, e2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First element selection:\n" << e1 << "\n"
                << "   Second element selection:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching element selections (different size)
      {
         auto e1 = blaze::elements( vec_, { 5UL, 3UL, 1UL } );
         auto e2 = blaze::elements( vec_, { 5UL, 3UL } );

         if( blaze::isSame( e1, e2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First element selection:\n" << e1 << "\n"
                << "   Second element selection:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching element selections (different order)
      {
         auto e1 = blaze::elements( vec_, { 5UL, 3UL, 1UL } );
         auto e2 = blaze::elements( vec_, { 5UL, 1UL, 3UL } );

         if( blaze::isSame( e1, e2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First element selection:\n" << e1 << "\n"
                << "   Second element selection:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-based tests
   //=====================================================================================

   {
      test_ = "isSame() function (row-based)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 2, 3 },
                                                           { 4, 5, 6 },
                                                           { 7, 8, 9 } };

      // isSame with row and matching element selection
      {
         auto r = blaze::row( mat, 1UL );
         auto e = blaze::elements( r, { 0UL, 1UL, 2UL } );

         if( blaze::isSame( e, r ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( r, e ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with row and non-matching element selection (different size)
      {
         auto r = blaze::row( mat, 1UL );
         auto e = blaze::elements( r, { 0UL, 1UL } );

         if( blaze::isSame( e, r ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( r, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with row and non-matching element selection (different order)
      {
         auto r = blaze::row( mat, 1UL );
         auto e = blaze::elements( r, { 0UL, 2UL, 1UL } );

         if( blaze::isSame( e, r ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( r, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Row:\n" << r << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with subrow and matching element selection
      {
         auto r = blaze::row( mat, 1UL );
         auto e = blaze::elements( r, { 1UL, 2UL } );
         auto s = blaze::subvector( r, 1UL, 2UL );

         if( blaze::isSame( e, s ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subrow:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( s, e ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subrow:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with subrow and non-matching element selection (different size)
      {
         auto r = blaze::row( mat, 1UL );
         auto e = blaze::elements( r, { 0UL, 1UL, 2UL } );
         auto s = blaze::subvector( r, 1UL, 2UL );

         if( blaze::isSame( e, s ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subrow:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( s, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subrow:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with subrow and non-matching element selection (different order)
      {
         auto r = blaze::row( mat, 1UL );
         auto e = blaze::elements( r, { 2UL, 1UL } );
         auto s = blaze::subvector( r, 1UL, 2UL );

         if( blaze::isSame( e, s ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subrow:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( s, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subrow:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching element selections
      {
         auto r  = blaze::row( mat, 1UL );
         auto e1 = blaze::elements( r, { 1UL, 2UL } );
         auto e2 = blaze::elements( r, { 1UL, 2UL } );

         if( blaze::isSame( e1, e2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First element selection:\n" << e1 << "\n"
                << "   Second element selection:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching element selections (different size)
      {
         auto r  = blaze::row( mat, 1UL );
         auto e1 = blaze::elements( r, { 1UL, 2UL } );
         auto e2 = blaze::elements( r, { 0UL, 1UL, 2UL } );

         if( blaze::isSame( e1, e2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First element selection:\n" << e1 << "\n"
                << "   Second element selection:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching element selections (different order)
      {
         auto r  = blaze::row( mat, 1UL );
         auto e1 = blaze::elements( r, { 1UL, 2UL } );
         auto e2 = blaze::elements( r, { 2UL, 1UL } );

         if( blaze::isSame( e1, e2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First element selection:\n" << e1 << "\n"
                << "   Second element selection:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-based tests
   //=====================================================================================

   {
      test_ = "isSame() function (column-based)";

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 2, 3 },
                                                           { 4, 5, 6 },
                                                           { 7, 8, 9 } };

      // isSame with column and matching element selection
      {
         auto c = blaze::column( mat, 1UL );
         auto e = blaze::elements( c, { 0UL, 1UL, 2UL } );

         if( blaze::isSame( e, c ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( c, e ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and non-matching element selection (different size)
      {
         auto c = blaze::column( mat, 1UL );
         auto e = blaze::elements( c, { 0UL, 1UL } );

         if( blaze::isSame( e, c ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( c, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with column and non-matching element selection (different order)
      {
         auto c = blaze::column( mat, 1UL );
         auto e = blaze::elements( c, { 0UL, 2UL, 1UL } );

         if( blaze::isSame( e, c ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( c, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Column:\n" << c << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with subcolumn and matching element selection
      {
         auto c = blaze::column( mat, 1UL );
         auto e = blaze::elements( c, { 1UL, 2UL } );
         auto s = blaze::subvector( c, 1UL, 2UL );

         if( blaze::isSame( e, s ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subcolumn:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( s, e ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subcolumn:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with subcolumn and non-matching element selection (different size)
      {
         auto c = blaze::column( mat, 1UL );
         auto e = blaze::elements( c, { 0UL, 1UL, 2UL } );
         auto s = blaze::subvector( c, 1UL, 2UL );

         if( blaze::isSame( e, s ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subcolumn:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( s, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subcolumn:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with subcolumn and non-matching element selection (different order)
      {
         auto c = blaze::column( mat, 1UL );
         auto e = blaze::elements( c, { 2UL, 1UL } );
         auto s = blaze::subvector( c, 1UL, 2UL );

         if( blaze::isSame( e, s ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subcolumn:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( blaze::isSame( s, e ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   Subcolumn:\n" << s << "\n"
                << "   Element selection:\n" << e << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with matching element selections
      {
         auto c  = blaze::column( mat, 1UL );
         auto e1 = blaze::elements( c, { 1UL, 2UL } );
         auto e2 = blaze::elements( c, { 1UL, 2UL } );

         if( blaze::isSame( e1, e2 ) == false ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First element selection:\n" << e1 << "\n"
                << "   Second element selection:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching element selections (different size)
      {
         auto c  = blaze::column( mat, 1UL );
         auto e1 = blaze::elements( c, { 1UL, 2UL } );
         auto e2 = blaze::elements( c, { 0UL, 1UL, 2UL } );

         if( blaze::isSame( e1, e2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First element selection:\n" << e1 << "\n"
                << "   Second element selection:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // isSame with non-matching element selections (different order)
      {
         auto c  = blaze::column( mat, 1UL );
         auto e1 = blaze::elements( c, { 1UL, 2UL } );
         auto e2 = blaze::elements( c, { 2UL, 1UL } );

         if( blaze::isSame( e1, e2 ) == true ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid isSame evaluation\n"
                << " Details:\n"
                << "   First element selection:\n" << e1 << "\n"
                << "   Second element selection:\n" << e2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c subvector() function with the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c subvector() function used with the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testSubvector()
{
   test_ = "subvector() function";

   initialize();

   {
      auto e = blaze::elements( vec_, { 1UL, 3UL, 5UL, 2UL, 4UL, 6UL } );
      auto s = blaze::subvector( e, 1UL, 4UL );

      if( s[0] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subscript operator access failed\n"
             << " Details:\n"
             << "   Result: " << s[0] << "\n"
             << "   Expected result: -2\n";
         throw std::runtime_error( oss.str() );
      }

      if( s.begin()->value() != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator access failed\n"
             << " Details:\n"
             << "   Result: " << s.begin()->value() << "\n"
             << "   Expected result: -2\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      auto e = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL, 5UL, 6UL } );
      auto s = blaze::subvector( e, 6UL, 4UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << s << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}

   try {
      auto e = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL, 5UL, 6UL } );
      auto s = blaze::subvector( e, 2UL, 5UL );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Setup of out-of-bounds subvector succeeded\n"
          << " Details:\n"
          << "   Result:\n" << s << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c elements() function with the Elements class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c elements() function used with the Elements
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseTest::testElements()
{
   //=====================================================================================
   // Setup via index_sequence
   //=====================================================================================

   {
      test_ = "elements() function (index_sequence)";

      initialize();

      {
         auto e1 = blaze::elements( vec_, { 1UL, 3UL, 5UL, 2UL, 4UL, 6UL } );
         auto e2 = blaze::elements( e1  , { 1UL, 2UL, 3UL, 4UL } );

         if( e2[0] != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e2[0] << "\n"
                << "   Expected result: -2\n";
            throw std::runtime_error( oss.str() );
         }

         if( e2.begin()->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << e2.begin()->value() << "\n"
                << "   Expected result: -2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         auto e1 = blaze::elements( vec_, { 3UL, 6UL } );
         auto e2 = blaze::elements( e1  , { 1UL, 1UL, 1UL } );

         if( e2[0] != 4 || e2[1] != 4 || e2[2] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e2[0] << "\n"
                << "   Expected result: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( e2.begin()->value() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << e2.begin()->value() << "\n"
                << "   Expected result: 4\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         auto e1 = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL, 5UL, 6UL } );
         auto e2 = blaze::elements( e1  , { 6UL } );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds elements succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e2 << "\n";
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
         std::array<int,4UL> indices{ 1UL, 2UL, 3UL, 4UL };

         auto e1 = blaze::elements( vec_, { 1UL, 3UL, 5UL, 2UL, 4UL, 6UL } );
         auto e2 = blaze::elements( e1  , indices );

         if( e2[0] != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e2[0] << "\n"
                << "   Expected result: -2\n";
            throw std::runtime_error( oss.str() );
         }

         if( e2.begin()->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << e2.begin()->value() << "\n"
                << "   Expected result: -2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         std::array<int,4UL> indices{ 1UL, 1UL, 1UL };

         auto e1 = blaze::elements( vec_, { 3UL, 6UL } );
         auto e2 = blaze::elements( e1  , indices );

         if( e2[0] != 4 || e2[1] != 4 || e2[2] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e2[0] << "\n"
                << "   Expected result: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( e2.begin()->value() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << e2.begin()->value() << "\n"
                << "   Expected result: 4\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         std::array<int,1UL> indices{ 6UL };

         auto e1 = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL, 5UL, 6UL } );
         auto e2 = blaze::elements( e1  , indices );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds elements succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e2 << "\n";
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
         auto e1 = blaze::elements( vec_, { 1UL, 3UL, 5UL, 2UL, 4UL, 6UL } );
         auto e2 = blaze::elements( e1  , []( size_t i ){ return i+1UL; }, 4UL );

         if( e2[0] != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e2[0] << "\n"
                << "   Expected result: -2\n";
            throw std::runtime_error( oss.str() );
         }

         if( e2.begin()->value() != -2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << e2.begin()->value() << "\n"
                << "   Expected result: -2\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         auto e1 = blaze::elements( vec_, { 3UL, 6UL } );
         auto e2 = blaze::elements( e1  , []( size_t ){ return 1UL; }, 3UL );

         if( e2[0] != 4 || e2[1] != 4 || e2[2] != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Subscript operator access failed\n"
                << " Details:\n"
                << "   Result: " << e2[0] << "\n"
                << "   Expected result: 4\n";
            throw std::runtime_error( oss.str() );
         }

         if( e2.begin()->value() != 4 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Iterator access failed\n"
                << " Details:\n"
                << "   Result: " << e2.begin()->value() << "\n"
                << "   Expected result: 4\n";
            throw std::runtime_error( oss.str() );
         }
      }

      try {
         auto e1 = blaze::elements( vec_, { 1UL, 2UL, 3UL, 4UL, 5UL, 6UL } );
         auto e2 = blaze::elements( e1  , []( size_t ){ return 6UL; }, 1UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds elements succeeded\n"
             << " Details:\n"
             << "   Result:\n" << e2 << "\n";
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
void SparseTest::initialize()
{
   // Initializing the compressed row vector
   vec_.reset();
   vec_[1] =  1;
   vec_[3] = -2;
   vec_[4] = -3;
   vec_[6] =  4;
}
//*************************************************************************************************

} // namespace elements

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
   std::cout << "   Running Elements sparse test..." << std::endl;

   try
   {
      RUN_ELEMENTS_SPARSE_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during Elements sparse test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
