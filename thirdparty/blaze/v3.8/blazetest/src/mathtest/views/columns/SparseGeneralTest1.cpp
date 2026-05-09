//=================================================================================================
/*!
//  \file src/mathtest/views/columns/SparseGeneralTest1.cpp
//  \brief Source file for the Columns sparse general test (part 1)
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
#include <blaze/math/DynamicVector.h>
#include <blazetest/mathtest/views/columns/SparseGeneralTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace views {

namespace columns {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the Columns sparse general test.
//
// \exception std::runtime_error Operation error detected.
*/
SparseGeneralTest::SparseGeneralTest()
   : mat_ ( 4UL, 5UL )
   , tmat_( 4UL, 5UL )
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
/*!\brief Test of the Columns constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the Columns specialization. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testConstructors()
{
   using blaze::index_sequence;
   using blaze::initializer_list;


   //=====================================================================================
   // Row-major setup via index_sequence
   //=====================================================================================

   {
      test_ = "Row-major Columns constructor (index_sequence)";

      initialize();

      // Setup of a regular column selection
      {
         auto cs = blaze::columns( mat_, index_sequence<0,4,2>() );

         if( cs.rows() != mat_.rows() || cs.columns() != 3UL ||
             cs(0,0) != mat_(0,0) || cs(0,1) != mat_(0,4) || cs(0,2) != mat_(0,2) ||
             cs(1,0) != mat_(1,0) || cs(1,1) != mat_(1,4) || cs(1,2) != mat_(1,2) ||
             cs(2,0) != mat_(2,0) || cs(2,1) != mat_(2,4) || cs(2,2) != mat_(2,2) ||
             cs(3,0) != mat_(3,0) || cs(3,1) != mat_(3,4) || cs(3,2) != mat_(3,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds column selection
      try {
         auto cs = blaze::columns( mat_, index_sequence<5>() );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds column selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of a column selection on a compile-time column selection
      {
         auto cs1 = blaze::columns( mat_, index_sequence<0,4,2>() );
         auto cs2 = blaze::columns( cs1, index_sequence<2,1>() );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an explicit column selection
      {
         auto cs1 = blaze::columns( mat_, { 0, 4, 2 } );
         auto cs2 = blaze::columns( cs1, index_sequence<2,1>() );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an implicit column selection
      {
         const std::array<size_t,3UL> indices{ 0, 4, 2 };
         auto cs1 = blaze::columns( mat_, [indices]( size_t i ){ return indices[i]; }, 3UL );
         auto cs2 = blaze::columns( cs1, index_sequence<2,1>() );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major setup via initializer_list
   //=====================================================================================

   {
      test_ = "Row-major Columns constructor (initializer_list)";

      initialize();

      // Setup of empty column selection
      {
         std::initializer_list<size_t> indices{};
         auto cs = blaze::columns( mat_, indices );

         if( cs.rows() != mat_.rows() || cs.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of empty column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a regular column selection
      {
         auto cs = blaze::columns( mat_, { 0, 4, 2 } );

         if( cs.rows() != mat_.rows() || cs.columns() != 3UL ||
             cs(0,0) != mat_(0,0) || cs(0,1) != mat_(0,4) || cs(0,2) != mat_(0,2) ||
             cs(1,0) != mat_(1,0) || cs(1,1) != mat_(1,4) || cs(1,2) != mat_(1,2) ||
             cs(2,0) != mat_(2,0) || cs(2,1) != mat_(2,4) || cs(2,2) != mat_(2,2) ||
             cs(3,0) != mat_(3,0) || cs(3,1) != mat_(3,4) || cs(3,2) != mat_(3,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds column selection
      try {
         auto cs = blaze::columns( mat_, { 5 } );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds column selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of a column selection on a compile-time column selection
      {
         auto cs1 = blaze::columns( mat_, index_sequence<0,4,2>() );
         auto cs2 = blaze::columns( cs1, { 2, 1 } );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an explicit column selection
      {
         auto cs1 = blaze::columns( mat_, { 0, 4, 2 } );
         auto cs2 = blaze::columns( cs1, { 2, 1 } );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an implicit column selection
      {
         const std::array<size_t,3UL> indices{ 0, 4, 2 };
         auto cs1 = blaze::columns( mat_, [indices]( size_t i ){ return indices[i]; }, 3UL );
         auto cs2 = blaze::columns( cs1, { 2, 1 } );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major setup via std::vector
   //=====================================================================================

   {
      test_ = "Row-major Columns constructor (std::vector)";

      initialize();

      // Setup of empty column selection
      {
         std::vector<size_t> indices;
         auto cs = blaze::columns( mat_, indices );

         if( cs.rows() != mat_.rows() || cs.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of empty column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a regular column selection
      {
         const std::vector<size_t> indices{ 0, 4, 2 };
         auto cs = blaze::columns( mat_, indices );

         if( cs.rows() != mat_.rows() || cs.columns() != 3UL ||
             cs(0,0) != mat_(0,0) || cs(0,1) != mat_(0,4) || cs(0,2) != mat_(0,2) ||
             cs(1,0) != mat_(1,0) || cs(1,1) != mat_(1,4) || cs(1,2) != mat_(1,2) ||
             cs(2,0) != mat_(2,0) || cs(2,1) != mat_(2,4) || cs(2,2) != mat_(2,2) ||
             cs(3,0) != mat_(3,0) || cs(3,1) != mat_(3,4) || cs(3,2) != mat_(3,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds column selection
      try {
         std::vector<size_t> indices{ 5 };
         auto cs = blaze::columns( mat_, indices );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds column selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of a column selection on a compile-time column selection
      {
         auto cs1 = blaze::columns( mat_, index_sequence<0,4,2>() );

         const std::vector<size_t> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an explicit column selection
      {
         auto cs1 = blaze::columns( mat_, { 0, 4, 2 } );

         const std::vector<size_t> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an implicit column selection
      {
         const std::array<size_t,3UL> indices1{ 0, 4, 2 };
         auto cs1 = blaze::columns( mat_, [indices1]( size_t i ){ return indices1[i]; }, 3UL );

         const std::vector<size_t> indices2{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices2 );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major setup via std::array
   //=====================================================================================

   {
      test_ = "Row-major Columns constructor (std::array)";

      initialize();

      // Setup of a regular column selection
      {
         const std::array<size_t,3UL> indices{ 0, 4, 2 };
         auto cs = blaze::columns( mat_, indices );

         if( cs.rows() != mat_.rows() || cs.columns() != 3UL ||
             cs(0,0) != mat_(0,0) || cs(0,1) != mat_(0,4) || cs(0,2) != mat_(0,2) ||
             cs(1,0) != mat_(1,0) || cs(1,1) != mat_(1,4) || cs(1,2) != mat_(1,2) ||
             cs(2,0) != mat_(2,0) || cs(2,1) != mat_(2,4) || cs(2,2) != mat_(2,2) ||
             cs(3,0) != mat_(3,0) || cs(3,1) != mat_(3,4) || cs(3,2) != mat_(3,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds column selection
      try {
         std::array<size_t,1UL> indices{ 5 };
         auto cs = blaze::columns( mat_, indices );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds column selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of a column selection on a compile-time column selection
      {
         auto cs1 = blaze::columns( mat_, index_sequence<0,4,2>() );

         const std::array<size_t,2UL> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an explicit column selection
      {
         auto cs1 = blaze::columns( mat_, { 0, 4, 2 } );

         const std::array<size_t,2UL> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an implicit column selection
      {
         const std::array<size_t,3UL> indices1{ 0, 4, 2 };
         auto cs1 = blaze::columns( mat_, [indices1]( size_t i ){ return indices1[i]; }, 3UL );

         const std::array<size_t,2UL> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major setup via lambda expression
   //=====================================================================================

   {
      test_ = "Row-major Columns constructor (lambda expression)";

      initialize();

      // Setup of empty column selection
      {
         auto cs = blaze::columns( mat_, []( size_t ){ return 0UL; }, 0UL );

         if( cs.rows() != mat_.rows() || cs.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of empty column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a regular column selection
      {
         const std::array<size_t,3UL> indices{ 0, 4, 2 };
         auto cs = blaze::columns( mat_, [indices]( size_t i ){ return indices[i]; }, 3UL );

         if( cs.rows() != mat_.rows() || cs.columns() != 3UL ||
             cs(0,0) != mat_(0,0) || cs(0,1) != mat_(0,4) || cs(0,2) != mat_(0,2) ||
             cs(1,0) != mat_(1,0) || cs(1,1) != mat_(1,4) || cs(1,2) != mat_(1,2) ||
             cs(2,0) != mat_(2,0) || cs(2,1) != mat_(2,4) || cs(2,2) != mat_(2,2) ||
             cs(3,0) != mat_(3,0) || cs(3,1) != mat_(3,4) || cs(3,2) != mat_(3,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds column selection
      try {
         auto cs = blaze::columns( mat_, []( size_t ){ return 5UL; }, 1UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds column selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of a column selection on a compile-time column selection
      {
         auto cs1 = blaze::columns( mat_, index_sequence<0,4,2>() );

         const std::array<size_t,2UL> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, [indices]( size_t i ){ return indices[i]; }, 2UL );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an explicit column selection
      {
         auto cs1 = blaze::columns( mat_, { 0, 4, 2 } );

         const std::array<size_t,2UL> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, [indices]( size_t i ){ return indices[i]; }, 2UL );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an implicit column selection
      {
         const std::array<size_t,3UL> indices1{ 0, 4, 2 };
         auto cs1 = blaze::columns( mat_, [indices1]( size_t i ){ return indices1[i]; }, 3UL );

         const std::array<size_t,2UL> indices2{ 2, 1 };
         auto cs2 = blaze::columns( cs1, [indices2]( size_t i ){ return indices2[i]; }, 2UL );

         if( cs2.rows() != mat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != mat_(0,2) || cs2(0,1) != mat_(0,4) ||
             cs2(1,0) != mat_(1,2) || cs2(1,1) != mat_(1,4) ||
             cs2(2,0) != mat_(2,2) || cs2(2,1) != mat_(2,4) ||
             cs2(3,0) != mat_(3,2) || cs2(3,1) != mat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major setup of random in-bounds element selection
   //=====================================================================================

   {
      test_ = "Column-major Columns constructor (stress test)";

      initialize();

      for( size_t rep=0UL; rep<100UL; ++rep )
      {
         blaze::DynamicVector<size_t> indices( blaze::rand<size_t>( 1UL, 20UL ) );
         randomize( indices, 0UL, mat_.rows()-1UL );
         auto cs = blaze::columns( mat_, indices.data(), indices.size() );

         for( size_t i=0UL; i<cs.rows(); ++i ) {
            for( size_t j=0UL; j<cs.columns(); ++j ) {
               if( cs(i,j) != mat_(i,indices[j]) ) {
                  std::ostringstream oss;
                  oss << " Test: " << test_ << "\n"
                      << " Error: Setup of column selection failed\n"
                      << " Details:\n"
                      << "   Indices:\n" << indices << "\n"
                      << "   Column selection:\n" << cs << "\n"
                      << "   Matrix:\n" << mat_ << "\n";
                  throw std::runtime_error( oss.str() );
               }
            }
         }
      }
   }


   //=====================================================================================
   // Column-major setup via index_sequence
   //=====================================================================================

   {
      test_ = "Column-major Columns constructor (index_sequence)";

      initialize();

      // Setup of a regular column selection
      {
         auto cs = blaze::columns( tmat_, index_sequence<0,4,2>() );

         if( cs.rows() != tmat_.rows() || cs.columns() != 3UL ||
             cs(0,0) != tmat_(0,0) || cs(0,1) != tmat_(0,4) || cs(0,2) != tmat_(0,2) ||
             cs(1,0) != tmat_(1,0) || cs(1,1) != tmat_(1,4) || cs(1,2) != tmat_(1,2) ||
             cs(2,0) != tmat_(2,0) || cs(2,1) != tmat_(2,4) || cs(2,2) != tmat_(2,2) ||
             cs(3,0) != tmat_(3,0) || cs(3,1) != tmat_(3,4) || cs(3,2) != tmat_(3,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds column selection
      try {
         auto cs = blaze::columns( tmat_, index_sequence<5>() );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds column selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of a column selection on a compile-time column selection
      {
         auto cs1 = blaze::columns( tmat_, index_sequence<0,4,2>() );
         auto cs2 = blaze::columns( cs1, index_sequence<2,1>() );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an explicit column selection
      {
         auto cs1 = blaze::columns( tmat_, { 0, 4, 2 } );
         auto cs2 = blaze::columns( cs1, index_sequence<2,1>() );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an implicit column selection
      {
         const std::array<size_t,3UL> indices{ 0, 4, 2 };
         auto cs1 = blaze::columns( tmat_, [indices]( size_t i ){ return indices[i]; }, 3UL );
         auto cs2 = blaze::columns( cs1, index_sequence<2,1>() );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major setup via initializer_list
   //=====================================================================================

   {
      test_ = "Column-major Columns constructor (initializer_list)";

      initialize();

      // Setup of empty column selection
      {
         std::initializer_list<size_t> indices{};
         auto cs = blaze::columns( tmat_, indices );

         if( cs.rows() != tmat_.rows() || cs.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of empty column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a regular column selection
      {
         auto cs = blaze::columns( tmat_, { 0, 4, 2 } );

         if( cs.rows() != tmat_.rows() || cs.columns() != 3UL ||
             cs(0,0) != tmat_(0,0) || cs(0,1) != tmat_(0,4) || cs(0,2) != tmat_(0,2) ||
             cs(1,0) != tmat_(1,0) || cs(1,1) != tmat_(1,4) || cs(1,2) != tmat_(1,2) ||
             cs(2,0) != tmat_(2,0) || cs(2,1) != tmat_(2,4) || cs(2,2) != tmat_(2,2) ||
             cs(3,0) != tmat_(3,0) || cs(3,1) != tmat_(3,4) || cs(3,2) != tmat_(3,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds column selection
      try {
         auto cs = blaze::columns( tmat_, { 5 } );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds column selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of a column selection on a compile-time column selection
      {
         auto cs1 = blaze::columns( tmat_, index_sequence<0,4,2>() );
         auto cs2 = blaze::columns( cs1, { 2, 1 } );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an explicit column selection
      {
         auto cs1 = blaze::columns( tmat_, { 0, 4, 2 } );
         auto cs2 = blaze::columns( cs1, { 2, 1 } );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an implicit column selection
      {
         const std::array<size_t,3UL> indices{ 0, 4, 2 };
         auto cs1 = blaze::columns( tmat_, [indices]( size_t i ){ return indices[i]; }, 3UL );
         auto cs2 = blaze::columns( cs1, { 2, 1 } );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major setup via std::vector
   //=====================================================================================

   {
      test_ = "Column-major Columns constructor (std::vector)";

      initialize();

      // Setup of empty column selection
      {
         std::vector<size_t> indices;
         auto cs = blaze::columns( tmat_, indices );

         if( cs.rows() != tmat_.rows() || cs.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of empty column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a regular column selection
      {
         const std::vector<size_t> indices{ 0, 4, 2 };
         auto cs = blaze::columns( tmat_, indices );

         if( cs.rows() != tmat_.rows() || cs.columns() != 3UL ||
             cs(0,0) != tmat_(0,0) || cs(0,1) != tmat_(0,4) || cs(0,2) != tmat_(0,2) ||
             cs(1,0) != tmat_(1,0) || cs(1,1) != tmat_(1,4) || cs(1,2) != tmat_(1,2) ||
             cs(2,0) != tmat_(2,0) || cs(2,1) != tmat_(2,4) || cs(2,2) != tmat_(2,2) ||
             cs(3,0) != tmat_(3,0) || cs(3,1) != tmat_(3,4) || cs(3,2) != tmat_(3,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds column selection
      try {
         std::vector<size_t> indices{ 5 };
         auto cs = blaze::columns( tmat_, indices );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds column selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of a column selection on a compile-time column selection
      {
         auto cs1 = blaze::columns( tmat_, index_sequence<0,4,2>() );

         const std::vector<size_t> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an explicit column selection
      {
         auto cs1 = blaze::columns( tmat_, { 0, 4, 2 } );

         const std::vector<size_t> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an implicit column selection
      {
         const std::array<size_t,3UL> indices1{ 0, 4, 2 };
         auto cs1 = blaze::columns( tmat_, [indices1]( size_t i ){ return indices1[i]; }, 3UL );

         const std::vector<size_t> indices2{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices2 );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major setup via std::array
   //=====================================================================================

   {
      test_ = "Column-major Columns constructor (std::array)";

      initialize();

      // Setup of a regular column selection
      {
         const std::array<size_t,3UL> indices{ 0, 4, 2 };
         auto cs = blaze::columns( tmat_, indices );

         if( cs.rows() != tmat_.rows() || cs.columns() != 3UL ||
             cs(0,0) != tmat_(0,0) || cs(0,1) != tmat_(0,4) || cs(0,2) != tmat_(0,2) ||
             cs(1,0) != tmat_(1,0) || cs(1,1) != tmat_(1,4) || cs(1,2) != tmat_(1,2) ||
             cs(2,0) != tmat_(2,0) || cs(2,1) != tmat_(2,4) || cs(2,2) != tmat_(2,2) ||
             cs(3,0) != tmat_(3,0) || cs(3,1) != tmat_(3,4) || cs(3,2) != tmat_(3,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds column selection
      try {
         std::array<size_t,1UL> indices{ 5 };
         auto cs = blaze::columns( tmat_, indices );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds column selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of a column selection on a compile-time column selection
      {
         auto cs1 = blaze::columns( tmat_, index_sequence<0,4,2>() );

         const std::array<size_t,2UL> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an explicit column selection
      {
         auto cs1 = blaze::columns( tmat_, { 0, 4, 2 } );

         const std::array<size_t,2UL> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an implicit column selection
      {
         const std::array<size_t,3UL> indices1{ 0, 4, 2 };
         auto cs1 = blaze::columns( tmat_, [indices1]( size_t i ){ return indices1[i]; }, 3UL );

         const std::array<size_t,2UL> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, indices );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major setup via lambda expression
   //=====================================================================================

   {
      test_ = "Column-major Columns constructor (lambda expression)";

      initialize();

      // Setup of empty column selection
      {
         auto cs = blaze::columns( tmat_, []( size_t ){ return 0UL; }, 0UL );

         if( cs.rows() != tmat_.rows() || cs.columns() != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of empty column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a regular column selection
      {
         const std::array<size_t,3UL> indices{ 0, 4, 2 };
         auto cs = blaze::columns( tmat_, [indices]( size_t i ){ return indices[i]; }, 3UL );

         if( cs.rows() != tmat_.rows() || cs.columns() != 3UL ||
             cs(0,0) != tmat_(0,0) || cs(0,1) != tmat_(0,4) || cs(0,2) != tmat_(0,2) ||
             cs(1,0) != tmat_(1,0) || cs(1,1) != tmat_(1,4) || cs(1,2) != tmat_(1,2) ||
             cs(2,0) != tmat_(2,0) || cs(2,1) != tmat_(2,4) || cs(2,2) != tmat_(2,2) ||
             cs(3,0) != tmat_(3,0) || cs(3,1) != tmat_(3,4) || cs(3,2) != tmat_(3,2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Trying to setup an out-of-bounds column selection
      try {
         auto cs = blaze::columns( tmat_, []( size_t ){ return 5UL; }, 1UL );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Setup of out-of-bounds column selection succeeded\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}

      // Setup of a column selection on a compile-time column selection
      {
         auto cs1 = blaze::columns( tmat_, index_sequence<0,4,2>() );

         const std::array<size_t,2UL> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, [indices]( size_t i ){ return indices[i]; }, 2UL );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an explicit column selection
      {
         auto cs1 = blaze::columns( tmat_, { 0, 4, 2 } );

         const std::array<size_t,2UL> indices{ 2, 1 };
         auto cs2 = blaze::columns( cs1, [indices]( size_t i ){ return indices[i]; }, 2UL );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      // Setup of a column selection on an implicit column selection
      {
         const std::array<size_t,3UL> indices1{ 0, 4, 2 };
         auto cs1 = blaze::columns( tmat_, [indices1]( size_t i ){ return indices1[i]; }, 3UL );

         const std::array<size_t,2UL> indices2{ 2, 1 };
         auto cs2 = blaze::columns( cs1, [indices2]( size_t i ){ return indices2[i]; }, 2UL );

         if( cs2.rows() != tmat_.rows() || cs2.columns() != 2UL ||
             cs2(0,0) != tmat_(0,2) || cs2(0,1) != tmat_(0,4) ||
             cs2(1,0) != tmat_(1,2) || cs2(1,1) != tmat_(1,4) ||
             cs2(2,0) != tmat_(2,2) || cs2(2,1) != tmat_(2,4) ||
             cs2(3,0) != tmat_(3,2) || cs2(3,1) != tmat_(3,4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Setup of column selection failed\n"
                << " Details:\n"
                << "   Result:\n" << cs2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major setup of random in-bounds element selection
   //=====================================================================================

   {
      test_ = "Column-major Columns constructor (stress test)";

      initialize();

      for( size_t rep=0UL; rep<100UL; ++rep )
      {
         blaze::DynamicVector<size_t> indices( blaze::rand<size_t>( 1UL, 20UL ) );
         randomize( indices, 0UL, tmat_.rows()-1UL );
         auto cs = blaze::columns( tmat_, indices.data(), indices.size() );

         for( size_t i=0UL; i<cs.rows(); ++i ) {
            for( size_t j=0UL; j<cs.columns(); ++j ) {
               if( cs(i,j) != tmat_(i,indices[j]) ) {
                  std::ostringstream oss;
                  oss << " Test: " << test_ << "\n"
                      << " Error: Setup of column selection failed\n"
                      << " Details:\n"
                      << "   Indices:\n" << indices << "\n"
                      << "   Column selection:\n" << cs << "\n"
                      << "   Matrix:\n" << tmat_ << "\n";
                  throw std::runtime_error( oss.str() );
               }
            }
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Columns assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the Columns specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testAssignment()
{
   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   {
      test_ = "Row-major Columns list assignment (complete list)";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );
      cs = { { 11, 0 }, { 0, 13 }, { 0, 14 }, { 12, 0 } };

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != -2 || mat_(0,3) != 11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 13 || mat_(1,2) !=  0 || mat_(1,3) !=  0 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 14 || mat_(2,2) != -3 || mat_(2,3) !=  0 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Columns list assignment (incomplete list)";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );
      cs = { { 11 }, { 0, 13 }, { 0, 14 }, { 12 } };

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != -2 || mat_(0,3) != 11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 13 || mat_(1,2) !=  0 || mat_(1,3) !=  0 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 14 || mat_(2,2) != -3 || mat_(2,3) !=  0 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major Columns copy assignment (no aliasing)";

      initialize();

      MT mat{ { 0, 11, 0, 13, 0 },
              { 0,  0, 0, 14, 0 },
              { 0, 12, 0, 15, 0 },
              { 0,  0, 0, 16, 0 } };

      auto cs = blaze::columns( mat, { 3UL, 1UL } );
      cs = blaze::columns( mat_, { 3UL, 1UL } );

      checkRows    ( cs , 4UL );
      checkColumns ( cs , 2UL );
      checkNonZeros( cs , 4UL );
      checkRows    ( mat, 4UL );
      checkColumns ( mat, 5UL );
      checkNonZeros( mat, 4UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) !=  4 || cs(1,1) != 1 ||
          cs(2,0) !=  5 || cs(2,1) != 0 ||
          cs(3,0) != -6 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n(  4  1 )\n(  5  0 )\n( -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 || mat(0,3) !=  0 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 1 || mat(1,2) != 0 || mat(1,3) !=  4 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) !=  5 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) != 0 || mat(3,2) != 0 || mat(3,3) != -6 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0  0  0 )\n"
                                     "( 0  1  0  4  0 )\n"
                                     "( 0  0  0  5  0 )\n"
                                     "( 0  0  0 -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Columns copy assignment (aliasing)";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 4UL } );
      cs = blaze::columns( mat_, { 2UL, 3UL } );

      checkRows    ( cs  , 4UL );
      checkColumns ( cs  , 2UL );
      checkNonZeros( cs  , 5UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 8UL );

      if( cs(0,0) != -2 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) !=  4 ||
          cs(2,0) != -3 || cs(2,1) !=  5 ||
          cs(3,0) !=  0 || cs(3,1) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -2  0 )\n(  0  4 )\n( -3  5 )\n(  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) != -2 || mat_(0,4) !=  0 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  0 || mat_(1,4) !=  4 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) != -3 || mat_(2,4) !=  5 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=  0 || mat_(3,4) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 -2  0 )\n"
                                     "( 0  1  0  0  4 )\n"
                                     "( 0  0 -3 -3  5 )\n"
                                     "( 0  0  0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                         {  0, 13 },
                                                         {  0, 14 },
                                                         { 12,  0 } };

      cs = mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != -2 || mat_(0,3) != 11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 13 || mat_(1,2) !=  0 || mat_(1,3) !=  0 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 14 || mat_(2,2) != -3 || mat_(2,3) !=  0 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs = mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != -2 || mat_(0,3) != 11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 13 || mat_(1,2) !=  0 || mat_(1,3) !=  0 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 14 || mat_(2,2) != -3 || mat_(2,3) !=  0 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs = mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != -2 || mat_(0,3) != 11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 13 || mat_(1,2) !=  0 || mat_(1,3) !=  0 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 14 || mat_(2,2) != -3 || mat_(2,3) !=  0 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                                 {  0, 13 },
                                                                 {  0, 14 },
                                                                 { 12,  0 } };

      cs = mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != -2 || mat_(0,3) != 11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 13 || mat_(1,2) !=  0 || mat_(1,3) !=  0 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 14 || mat_(2,2) != -3 || mat_(2,3) !=  0 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) != 12 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   {
      test_ = "Column-major Columns list assignment (complete list)";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );
      cs = { { 11, 0 }, { 0, 13 }, { 0, 14 }, { 12, 0 } };

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  4UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) != 11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 13 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 14 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Columns list assignment (incomplete list)";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );
      cs = { { 11 }, { 0, 13 }, { 0, 14 }, { 12 } };

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) != 11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 13 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 14 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major Columns copy assignment (no aliasing)";

      initialize();

      OMT mat{ { 0, 11, 0, 13, 0 },
               { 0,  0, 0, 14, 0 },
               { 0, 12, 0, 15, 0 },
               { 0,  0, 0, 16, 0 } };

      auto cs = blaze::columns( mat, { 3UL, 1UL } );
      cs = blaze::columns( tmat_, { 3UL, 1UL } );

      checkRows    ( cs , 4UL );
      checkColumns ( cs , 2UL );
      checkNonZeros( cs , 4UL );
      checkRows    ( mat, 4UL );
      checkColumns ( mat, 5UL );
      checkNonZeros( mat, 4UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) !=  4 || cs(1,1) != 1 ||
          cs(2,0) !=  5 || cs(2,1) != 0 ||
          cs(3,0) != -6 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n(  4  1 )\n(  5  0 )\n( -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 || mat(0,3) !=  0 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 1 || mat(1,2) != 0 || mat(1,3) !=  4 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) !=  5 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) != 0 || mat(3,2) != 0 || mat(3,3) != -6 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0  0  0 )\n"
                                     "( 0  1  0  4  0 )\n"
                                     "( 0  0  0  5  0 )\n"
                                     "( 0  0  0 -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Columns copy assignment (aliasing)";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 4UL } );
      cs = blaze::columns( tmat_, { 2UL, 3UL } );

      checkRows    ( cs   , 4UL );
      checkColumns ( cs   , 2UL );
      checkNonZeros( cs   , 5UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 8UL );

      if( cs(0,0) != -2 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) !=  4 ||
          cs(2,0) != -3 || cs(2,1) !=  5 ||
          cs(3,0) !=  0 || cs(3,1) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -2  0 )\n(  0  4 )\n( -3  5 )\n(  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) != -2 || tmat_(0,4) !=  0 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 || tmat_(1,4) !=  4 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) != -3 || tmat_(2,4) !=  5 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) !=  0 || tmat_(3,4) != -6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 -2  0 )\n"
                                     "( 0  1  0  0  4 )\n"
                                     "( 0  0 -3 -3  5 )\n"
                                     "( 0  0  0  0 -6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                           {  0, 13 },
                                                           {  0, 14 },
                                                           { 12,  0 } };

      cs = mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) != 11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 13 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 14 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs = mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) != 11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 13 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 14 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs = mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) != 11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 13 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 14 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                                 {  0, 13 },
                                                                 {  0, 14 },
                                                                 { 12,  0 } };

      cs = mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  4UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 10UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  0 || cs(1,1) != 13 ||
          cs(2,0) !=  0 || cs(2,1) != 14 ||
          cs(3,0) != 12 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) != 11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 13 || tmat_(1,2) !=  0 || tmat_(1,3) !=  0 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 14 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 13  0  0 -8 )\n"
                                     "( 0 14 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Columns addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the Columns specialization.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testAddAssign()
{
   //=====================================================================================
   // Row-major Columns addition assignment
   //=====================================================================================

   {
      test_ = "Row-major Columns addition assignment (no aliasing)";

      initialize();

      MT mat{ { 0, 11, 0, 13, 0 },
              { 0,  0, 0, 14, 0 },
              { 0, 12, 0, 15, 0 },
              { 0,  0, 0, 16, 0 } };

      auto cs = blaze::columns( mat, { 3UL, 1UL } );
      cs += blaze::columns( mat_, { 3UL, 1UL } );

      checkRows    ( cs , 4UL );
      checkColumns ( cs , 2UL );
      checkNonZeros( cs , 7UL );
      checkRows    ( mat, 4UL );
      checkColumns ( mat, 5UL );
      checkNonZeros( mat, 7UL );

      if( cs(0,0) != 13 || cs(0,1) != 11 ||
          cs(1,0) != 18 || cs(1,1) !=  1 ||
          cs(2,0) != 20 || cs(2,1) != 12 ||
          cs(3,0) != 10 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 13 11 )\n( 18  1 )\n( 20 12 )\n( 10  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 11 || mat(0,2) != 0 || mat(0,3) != 13 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) !=  1 || mat(1,2) != 0 || mat(1,3) != 18 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 12 || mat(2,2) != 0 || mat(2,3) != 20 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) !=  0 || mat(3,2) != 0 || mat(3,3) != 10 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 11  0 13  0 )\n"
                                     "( 0  1  0 18  0 )\n"
                                     "( 0 12  0 20  0 )\n"
                                     "( 0  0  0 10  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Columns addition assignment (aliasing)";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 4UL } );
      cs += blaze::columns( mat_, { 2UL, 3UL } );

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  8UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 11UL );

      if( cs(0,0) != -2 || cs(0,1) !=  7 ||
          cs(1,0) !=  4 || cs(1,1) != -4 ||
          cs(2,0) !=  2 || cs(2,1) != 14 ||
          cs(3,0) != -6 || cs(3,1) !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -2  7 )\n(  4 -4 )\n(  2 14 )\n( -6  4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) != -2 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -4 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  2 || mat_(2,4) != 14 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 -2  7 )\n"
                                     "( 0  1  0  4 -4 )\n"
                                     "( 0  0 -3  2 14 )\n"
                                     "( 0  0  0 -6  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix addition assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                           {  0, 13 },
                                                           {  0, 14 },
                                                           { 12,  0 } };

      cs += mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  6UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  4 || cs(1,1) != 14 ||
          cs(2,0) !=  5 || cs(2,1) != 14 ||
          cs(3,0) !=  6 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != -2 || mat_(0,3) != 11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 14 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) !=  6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 14  0  4 -8 )\n"
                                     "( 0 14 -3  5  9 )\n"
                                     "( 0  0  0  6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix addition assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs += mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  6UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  4 || cs(1,1) != 14 ||
          cs(2,0) !=  5 || cs(2,1) != 14 ||
          cs(3,0) !=  6 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != -2 || mat_(0,3) != 11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 14 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) !=  6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 14  0  4 -8 )\n"
                                     "( 0 14 -3  5  9 )\n"
                                     "( 0  0  0  6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix addition assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs += mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  6UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  4 || cs(1,1) != 14 ||
          cs(2,0) !=  5 || cs(2,1) != 14 ||
          cs(3,0) !=  6 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != -2 || mat_(0,3) != 11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 14 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) !=  6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 14  0  4 -8 )\n"
                                     "( 0 14 -3  5  9 )\n"
                                     "( 0  0  0  6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix addition assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                                 {  0, 13 },
                                                                 {  0, 14 },
                                                                 { 12,  0 } };

      cs += mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  6UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  4 || cs(1,1) != 14 ||
          cs(2,0) !=  5 || cs(2,1) != 14 ||
          cs(3,0) !=  6 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=  0 || mat_(0,2) != -2 || mat_(0,3) != 11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 14 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 14 || mat_(2,2) != -3 || mat_(2,3) !=  5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=  0 || mat_(3,2) !=  0 || mat_(3,3) !=  6 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 14  0  4 -8 )\n"
                                     "( 0 14 -3  5  9 )\n"
                                     "( 0  0  0  6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Columns addition assignment
   //=====================================================================================

   {
      test_ = "Column-major Columns addition assignment (no aliasing)";

      initialize();

      OMT mat{ { 0, 11, 0, 13, 0 },
               { 0,  0, 0, 14, 0 },
               { 0, 12, 0, 15, 0 },
               { 0,  0, 0, 16, 0 } };

      auto cs = blaze::columns( mat, { 3UL, 1UL } );
      cs += blaze::columns( tmat_, { 3UL, 1UL } );

      checkRows    ( cs , 4UL );
      checkColumns ( cs , 2UL );
      checkNonZeros( cs , 7UL );
      checkRows    ( mat, 4UL );
      checkColumns ( mat, 5UL );
      checkNonZeros( mat, 7UL );

      if( cs(0,0) != 13 || cs(0,1) != 11 ||
          cs(1,0) != 18 || cs(1,1) !=  1 ||
          cs(2,0) != 20 || cs(2,1) != 12 ||
          cs(3,0) != 10 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 13 11 )\n( 18  1 )\n( 20 12 )\n( 10  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 11 || mat(0,2) != 0 || mat(0,3) != 13 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) !=  1 || mat(1,2) != 0 || mat(1,3) != 18 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 12 || mat(2,2) != 0 || mat(2,3) != 20 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) !=  0 || mat(3,2) != 0 || mat(3,3) != 10 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 11  0 13  0 )\n"
                                     "( 0  1  0 18  0 )\n"
                                     "( 0 12  0 20  0 )\n"
                                     "( 0  0  0 10  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Columns addition assignment (aliasing)";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 4UL } );
      cs += blaze::columns( tmat_, { 2UL, 3UL } );

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  8UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 11UL );

      if( cs(0,0) != -2 || cs(0,1) !=  7 ||
          cs(1,0) !=  4 || cs(1,1) != -4 ||
          cs(2,0) !=  2 || cs(2,1) != 14 ||
          cs(3,0) != -6 || cs(3,1) !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -2  7 )\n(  4 -4 )\n(  2 14 )\n( -6  4 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) != -2 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -4 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  2 || tmat_(2,4) != 14 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) !=  4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 -2  7 )\n"
                                     "( 0  1  0  4 -4 )\n"
                                     "( 0  0 -3  2 14 )\n"
                                     "( 0  0  0 -6  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix addition assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                           {  0, 13 },
                                                           {  0, 14 },
                                                           { 12,  0 } };

      cs += mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  6UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  4 || cs(1,1) != 14 ||
          cs(2,0) !=  5 || cs(2,1) != 14 ||
          cs(3,0) !=  6 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) != 11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 14 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 14 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) !=  6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 14  0  4 -8 )\n"
                                     "( 0 14 -3  5  9 )\n"
                                     "( 0  0  0  6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix addition assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs += mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  6UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  4 || cs(1,1) != 14 ||
          cs(2,0) !=  5 || cs(2,1) != 14 ||
          cs(3,0) !=  6 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) != 11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 14 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 14 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) !=  6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 14  0  4 -8 )\n"
                                     "( 0 14 -3  5  9 )\n"
                                     "( 0  0  0  6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix addition assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs += mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  6UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  4 || cs(1,1) != 14 ||
          cs(2,0) !=  5 || cs(2,1) != 14 ||
          cs(3,0) !=  6 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) != 11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 14 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 14 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) !=  6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 14  0  4 -8 )\n"
                                     "( 0 14 -3  5  9 )\n"
                                     "( 0  0  0  6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix addition assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                                 {  0, 13 },
                                                                 {  0, 14 },
                                                                 { 12,  0 } };

      cs += mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  6UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( cs(0,0) != 11 || cs(0,1) !=  0 ||
          cs(1,0) !=  4 || cs(1,1) != 14 ||
          cs(2,0) !=  5 || cs(2,1) != 14 ||
          cs(3,0) !=  6 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 11  0 )\n(  0 13 )\n(  0 14 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=  0 || tmat_(0,2) != -2 || tmat_(0,3) != 11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 14 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 14 || tmat_(2,2) != -3 || tmat_(2,3) !=  5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=  0 || tmat_(3,2) !=  0 || tmat_(3,3) !=  6 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2 11  7 )\n"
                                     "( 0 14  0  4 -8 )\n"
                                     "( 0 14 -3  5  9 )\n"
                                     "( 0  0  0  6 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Columns subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the Columns
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testSubAssign()
{
   //=====================================================================================
   // Row-major Columns subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major Columns subtraction assignment (no aliasing)";

      initialize();

      MT mat{ { 0, 11, 0, 13, 0 },
              { 0,  0, 0, 14, 0 },
              { 0, 12, 0, 15, 0 },
              { 0,  0, 0, 16, 0 } };

      auto cs = blaze::columns( mat, { 3UL, 1UL } );
      cs -= blaze::columns( mat_, { 3UL, 1UL } );

      checkRows    ( cs , 4UL );
      checkColumns ( cs , 2UL );
      checkNonZeros( cs , 7UL );
      checkRows    ( mat, 4UL );
      checkColumns ( mat, 5UL );
      checkNonZeros( mat, 7UL );

      if( cs(0,0) != 13 || cs(0,1) != 11 ||
          cs(1,0) != 10 || cs(1,1) != -1 ||
          cs(2,0) != 10 || cs(2,1) != 12 ||
          cs(3,0) != 22 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 13 11 )\n( 10 -1 )\n( 10 12 )\n( 22  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 11 || mat(0,2) != 0 || mat(0,3) != 13 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) != -1 || mat(1,2) != 0 || mat(1,3) != 10 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 12 || mat(2,2) != 0 || mat(2,3) != 10 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) !=  0 || mat(3,2) != 0 || mat(3,3) != 22 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 11  0 13  0 )\n"
                                     "( 0 -1  0 10  0 )\n"
                                     "( 0 12  0 10  0 )\n"
                                     "( 0  0  0 22  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Columns subtraction assignment (aliasing)";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 4UL } );
      cs -= blaze::columns( mat_, { 2UL, 3UL } );

      checkRows    ( cs  , 4UL );
      checkColumns ( cs  , 2UL );
      checkNonZeros( cs  , 8UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 11UL );

      if( cs(0,0) !=  2 || cs(0,1) !=   7 ||
          cs(1,0) !=  4 || cs(1,1) != -12 ||
          cs(2,0) !=  8 || cs(2,1) !=   4 ||
          cs(3,0) != -6 || cs(3,1) !=  16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 13 11 )\n( 10 -1 )\n( 10 12 )\n( 22  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  2 || mat_(0,4) !=   7 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=  4 || mat_(1,4) != -12 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  8 || mat_(2,4) !=   4 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != -6 || mat_(3,4) !=  16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  2   7 )\n"
                                     "( 0  1  0  4 -12 )\n"
                                     "( 0  0 -3  8   4 )\n"
                                     "( 0  0  0 -6  16 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix subtraction assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                           {  0, 13 },
                                                           {  0, 14 },
                                                           { 12,  0 } };

      cs -= mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  6UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( cs(0,0) != -11 || cs(0,1) !=   0 ||
          cs(1,0) !=   4 || cs(1,1) != -12 ||
          cs(2,0) !=   5 || cs(2,1) != -14 ||
          cs(3,0) != -18 || cs(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -11   0 )\n(   4 -12 )\n(   5 -14 )\n( -18   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=   0 || mat_(0,2) != -2 || mat_(0,3) != -11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != -12 || mat_(1,2) !=  0 || mat_(1,3) !=   4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != -14 || mat_(2,2) != -3 || mat_(2,3) !=   5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=   0 || mat_(3,2) !=  0 || mat_(3,3) != -18 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0   0  -2 -11   7 )\n"
                                     "( 0 -12   0   4  -8 )\n"
                                     "( 0 -14  -3   5   9 )\n"
                                     "( 0   0   0 -18  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix subtraction assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs -= mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  6UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( cs(0,0) != -11 || cs(0,1) !=   0 ||
          cs(1,0) !=   4 || cs(1,1) != -12 ||
          cs(2,0) !=   5 || cs(2,1) != -14 ||
          cs(3,0) != -18 || cs(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -11   0 )\n(   4 -12 )\n(   5 -14 )\n( -18   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=   0 || mat_(0,2) != -2 || mat_(0,3) != -11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != -12 || mat_(1,2) !=  0 || mat_(1,3) !=   4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != -14 || mat_(2,2) != -3 || mat_(2,3) !=   5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=   0 || mat_(3,2) !=  0 || mat_(3,3) != -18 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0   0  -2 -11   7 )\n"
                                     "( 0 -12   0   4  -8 )\n"
                                     "( 0 -14  -3   5   9 )\n"
                                     "( 0   0   0 -18  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix subtraction assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs -= mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  6UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( cs(0,0) != -11 || cs(0,1) !=   0 ||
          cs(1,0) !=   4 || cs(1,1) != -12 ||
          cs(2,0) !=   5 || cs(2,1) != -14 ||
          cs(3,0) != -18 || cs(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -11   0 )\n(   4 -12 )\n(   5 -14 )\n( -18   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=   0 || mat_(0,2) != -2 || mat_(0,3) != -11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != -12 || mat_(1,2) !=  0 || mat_(1,3) !=   4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != -14 || mat_(2,2) != -3 || mat_(2,3) !=   5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=   0 || mat_(3,2) !=  0 || mat_(3,3) != -18 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0   0  -2 -11   7 )\n"
                                     "( 0 -12   0   4  -8 )\n"
                                     "( 0 -14  -3   5   9 )\n"
                                     "( 0   0   0 -18  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix subtraction assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                                 {  0, 13 },
                                                                 {  0, 14 },
                                                                 { 12,  0 } };

      cs -= mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  2UL );
      checkNonZeros( cs  ,  6UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( cs(0,0) != -11 || cs(0,1) !=   0 ||
          cs(1,0) !=   4 || cs(1,1) != -12 ||
          cs(2,0) !=   5 || cs(2,1) != -14 ||
          cs(3,0) != -18 || cs(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -11   0 )\n(   4 -12 )\n(   5 -14 )\n( -18   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) !=   0 || mat_(0,2) != -2 || mat_(0,3) != -11 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != -12 || mat_(1,2) !=  0 || mat_(1,3) !=   4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != -14 || mat_(2,2) != -3 || mat_(2,3) !=   5 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) !=   0 || mat_(3,2) !=  0 || mat_(3,3) != -18 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0   0  -2 -11   7 )\n"
                                     "( 0 -12   0   4  -8 )\n"
                                     "( 0 -14  -3   5   9 )\n"
                                     "( 0   0   0 -18  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Columns subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major Columns subtraction assignment (no aliasing)";

      initialize();

      OMT mat{ { 0, 11, 0, 13, 0 },
               { 0,  0, 0, 14, 0 },
               { 0, 12, 0, 15, 0 },
               { 0,  0, 0, 16, 0 } };

      auto cs = blaze::columns( mat, { 3UL, 1UL } );
      cs -= blaze::columns( mat_, { 3UL, 1UL } );

      checkRows    ( cs , 4UL );
      checkColumns ( cs , 2UL );
      checkNonZeros( cs , 7UL );
      checkRows    ( mat, 4UL );
      checkColumns ( mat, 5UL );
      checkNonZeros( mat, 7UL );

      if( cs(0,0) != 13 || cs(0,1) != 11 ||
          cs(1,0) != 10 || cs(1,1) != -1 ||
          cs(2,0) != 10 || cs(2,1) != 12 ||
          cs(3,0) != 22 || cs(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 13 11 )\n( 10 -1 )\n( 10 12 )\n( 22  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 11 || mat(0,2) != 0 || mat(0,3) != 13 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) != -1 || mat(1,2) != 0 || mat(1,3) != 10 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 12 || mat(2,2) != 0 || mat(2,3) != 10 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) !=  0 || mat(3,2) != 0 || mat(3,3) != 22 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 11  0 13  0 )\n"
                                     "( 0 -1  0 10  0 )\n"
                                     "( 0 12  0 10  0 )\n"
                                     "( 0  0  0 22  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Columns subtraction assignment (aliasing)";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 4UL } );
      cs -= blaze::columns( tmat_, { 2UL, 3UL } );

      checkRows    ( cs   , 4UL );
      checkColumns ( cs   , 2UL );
      checkNonZeros( cs   , 8UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 11UL );

      if( cs(0,0) !=  2 || cs(0,1) !=   7 ||
          cs(1,0) !=  4 || cs(1,1) != -12 ||
          cs(2,0) !=  8 || cs(2,1) !=   4 ||
          cs(3,0) != -6 || cs(3,1) !=  16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 13 11 )\n( 10 -1 )\n( 10 12 )\n( 22  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  2 || tmat_(0,4) !=   7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=  4 || tmat_(1,4) != -12 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  8 || tmat_(2,4) !=   4 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != -6 || tmat_(3,4) !=  16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  2   7 )\n"
                                     "( 0  1  0  4 -12 )\n"
                                     "( 0  0 -3  8   4 )\n"
                                     "( 0  0  0 -6  16 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix subtraction assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                           {  0, 13 },
                                                           {  0, 14 },
                                                           { 12,  0 } };

      cs -= mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  6UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( cs(0,0) != -11 || cs(0,1) !=   0 ||
          cs(1,0) !=   4 || cs(1,1) != -12 ||
          cs(2,0) !=   5 || cs(2,1) != -14 ||
          cs(3,0) != -18 || cs(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -11   0 )\n(   4 -12 )\n(   5 -14 )\n( -18   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=   0 || tmat_(0,2) != -2 || tmat_(0,3) != -11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -12 || tmat_(1,2) !=  0 || tmat_(1,3) !=   4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != -14 || tmat_(2,2) != -3 || tmat_(2,3) !=   5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=   0 || tmat_(3,2) !=  0 || tmat_(3,3) != -18 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0   0  -2 -11   7 )\n"
                                     "( 0 -12   0   4  -8 )\n"
                                     "( 0 -14  -3   5   9 )\n"
                                     "( 0   0   0 -18  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix subtraction assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs -= mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  6UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( cs(0,0) != -11 || cs(0,1) !=   0 ||
          cs(1,0) !=   4 || cs(1,1) != -12 ||
          cs(2,0) !=   5 || cs(2,1) != -14 ||
          cs(3,0) != -18 || cs(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -11   0 )\n(   4 -12 )\n(   5 -14 )\n( -18   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=   0 || tmat_(0,2) != -2 || tmat_(0,3) != -11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -12 || tmat_(1,2) !=  0 || tmat_(1,3) !=   4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != -14 || tmat_(2,2) != -3 || tmat_(2,3) !=   5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=   0 || tmat_(3,2) !=  0 || tmat_(3,3) != -18 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0   0  -2 -11   7 )\n"
                                     "( 0 -12   0   4  -8 )\n"
                                     "( 0 -14  -3   5   9 )\n"
                                     "( 0   0   0 -18  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix subtraction assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat{ { 11,  0 },
                                                              {  0, 13 },
                                                              {  0, 14 },
                                                              { 12,  0 } };

      cs -= mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  6UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( cs(0,0) != -11 || cs(0,1) !=   0 ||
          cs(1,0) !=   4 || cs(1,1) != -12 ||
          cs(2,0) !=   5 || cs(2,1) != -14 ||
          cs(3,0) != -18 || cs(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -11   0 )\n(   4 -12 )\n(   5 -14 )\n( -18   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=   0 || tmat_(0,2) != -2 || tmat_(0,3) != -11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -12 || tmat_(1,2) !=  0 || tmat_(1,3) !=   4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != -14 || tmat_(2,2) != -3 || tmat_(2,3) !=   5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=   0 || tmat_(3,2) !=  0 || tmat_(3,3) != -18 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0   0  -2 -11   7 )\n"
                                     "( 0 -12   0   4  -8 )\n"
                                     "( 0 -14  -3   5   9 )\n"
                                     "( 0   0   0 -18  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix subtraction assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat{ { 11,  0 },
                                                                 {  0, 13 },
                                                                 {  0, 14 },
                                                                 { 12,  0 } };

      cs -= mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  2UL );
      checkNonZeros( cs   ,  6UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( cs(0,0) != -11 || cs(0,1) !=   0 ||
          cs(1,0) !=   4 || cs(1,1) != -12 ||
          cs(2,0) !=   5 || cs(2,1) != -14 ||
          cs(3,0) != -18 || cs(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( -11   0 )\n(   4 -12 )\n(   5 -14 )\n( -18   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) !=   0 || tmat_(0,2) != -2 || tmat_(0,3) != -11 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != -12 || tmat_(1,2) !=  0 || tmat_(1,3) !=   4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != -14 || tmat_(2,2) != -3 || tmat_(2,3) !=   5 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) !=   0 || tmat_(3,2) !=  0 || tmat_(3,3) != -18 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0   0  -2 -11   7 )\n"
                                     "( 0 -12   0   4  -8 )\n"
                                     "( 0 -14  -3   5   9 )\n"
                                     "( 0   0   0 -18  10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Columns Schur product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment operators of the Columns
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testSchurAssign()
{
   //=====================================================================================
   // Row-major Columns Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major Columns Schur product assignment (no aliasing)";

      initialize();

      MT mat{ { 0, 1, 0, 4, 0 },
              { 0, 2, 0, 3, 0 },
              { 0, 3, 0, 2, 0 },
              { 0, 0, 0, 1, 0 } };

      auto cs = blaze::columns( mat, { 3UL, 1UL } );
      cs %= blaze::columns( mat_, { 3UL, 1UL } );

      checkRows    ( cs , 4UL );
      checkColumns ( cs , 2UL );
      checkNonZeros( cs , 4UL );
      checkRows    ( mat, 4UL );
      checkColumns ( mat, 5UL );
      checkNonZeros( mat, 4UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) != 12 || cs(1,1) != 2 ||
          cs(2,0) != 10 || cs(2,1) != 0 ||
          cs(3,0) != -6 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n( 12  2 )\n( 10  0 )\n( -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 || mat(0,3) !=  0 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 2 || mat(1,2) != 0 || mat(1,3) != 12 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) != 10 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) != 0 || mat(3,2) != 0 || mat(3,3) != -6 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0  0  0 )\n"
                                     "( 0  2  0 12  0 )\n"
                                     "( 0  0  0 10  0 )\n"
                                     "( 0  0  0 -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Columns Schur product assignment (aliasing)";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 4UL } );
      cs %= blaze::columns( mat_, { 2UL, 3UL } );

      checkRows    ( cs  , 4UL );
      checkColumns ( cs  , 2UL );
      checkNonZeros( cs  , 4UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 7UL );

      if( cs(0,0) !=   0 || cs(0,1) !=   0 ||
          cs(1,0) !=   0 || cs(1,1) != -32 ||
          cs(2,0) != -15 || cs(2,1) !=  45 ||
          cs(3,0) !=   0 || cs(3,1) != -60 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(   0   0 )\n(   0 -32 )\n( -15  45 )\n(   0 -60 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=   0 || mat_(0,4) !=   0 ||
          mat_(1,0) != 0 || mat_(1,1) != 1 || mat_(1,2) !=  0 || mat_(1,3) !=   0 || mat_(1,4) != -32 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) != -15 || mat_(2,4) !=  45 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) !=   0 || mat_(3,4) != -60 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0   0  -2   0   0 )\n"
                                     "( 0   1   0   0  32 )\n"
                                     "( 0   0  -3 -15  45 )\n"
                                     "( 0   0   0   0 -60 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix Schur product assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ {  0, 0 },
                                                           { -1, 2 },
                                                           {  0, 1 },
                                                           { -2, 0 } };

      cs %= mat;

      checkRows    ( cs  , 4UL );
      checkColumns ( cs  , 2UL );
      checkNonZeros( cs  , 3UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) != -4 || cs(1,1) != 2 ||
          cs(2,0) !=  0 || cs(2,1) != 0 ||
          cs(3,0) != 12 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n( -4  2 )\n(  0  0 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 2 || mat_(1,2) !=  0 || mat_(1,3) != -4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  0 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != 12 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  2  0 -4 -8 )\n"
                                     "( 0  0 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix Schur product assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ {  0, 0 },
                                                              { -1, 2 },
                                                              {  0, 1 },
                                                              { -2, 0 } };

      cs %= mat;

      checkRows    ( cs  , 4UL );
      checkColumns ( cs  , 2UL );
      checkNonZeros( cs  , 3UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) != -4 || cs(1,1) != 2 ||
          cs(2,0) !=  0 || cs(2,1) != 0 ||
          cs(3,0) != 12 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n( -4  2 )\n(  0  0 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 2 || mat_(1,2) !=  0 || mat_(1,3) != -4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  0 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != 12 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  2  0 -4 -8 )\n"
                                     "( 0  0 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix Schur product assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat{ {  0, 0 },
                                                              { -1, 2 },
                                                              {  0, 1 },
                                                              { -2, 0 } };

      cs %= mat;

      checkRows    ( cs  , 4UL );
      checkColumns ( cs  , 2UL );
      checkNonZeros( cs  , 3UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) != -4 || cs(1,1) != 2 ||
          cs(2,0) !=  0 || cs(2,1) != 0 ||
          cs(3,0) != 12 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n( -4  2 )\n(  0  0 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 2 || mat_(1,2) !=  0 || mat_(1,3) != -4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  0 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != 12 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  2  0 -4 -8 )\n"
                                     "( 0  0 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix Schur product assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat{ {  0, 0 },
                                                                 { -1, 2 },
                                                                 {  0, 1 },
                                                                 { -2, 0 } };

      cs %= mat;

      checkRows    ( cs  , 4UL );
      checkColumns ( cs  , 2UL );
      checkNonZeros( cs  , 3UL );
      checkRows    ( mat_, 4UL );
      checkColumns ( mat_, 5UL );
      checkNonZeros( mat_, 9UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) != -4 || cs(1,1) != 2 ||
          cs(2,0) !=  0 || cs(2,1) != 0 ||
          cs(3,0) != 12 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n( -4  2 )\n(  0  0 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != 0 || mat_(0,1) != 0 || mat_(0,2) != -2 || mat_(0,3) !=  0 || mat_(0,4) !=  7 ||
          mat_(1,0) != 0 || mat_(1,1) != 2 || mat_(1,2) !=  0 || mat_(1,3) != -4 || mat_(1,4) != -8 ||
          mat_(2,0) != 0 || mat_(2,1) != 0 || mat_(2,2) != -3 || mat_(2,3) !=  0 || mat_(2,4) !=  9 ||
          mat_(3,0) != 0 || mat_(3,1) != 0 || mat_(3,2) !=  0 || mat_(3,3) != 12 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  2  0 -4 -8 )\n"
                                     "( 0  0 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Columns Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major Columns Schur product assignment (no aliasing)";

      initialize();

      OMT mat{ { 0, 1, 0, 4, 0 },
               { 0, 2, 0, 3, 0 },
               { 0, 3, 0, 2, 0 },
               { 0, 0, 0, 1, 0 } };

      auto cs = blaze::columns( mat, { 3UL, 1UL } );
      cs %= blaze::columns( tmat_, { 3UL, 1UL } );

      checkRows    ( cs , 4UL );
      checkColumns ( cs , 2UL );
      checkNonZeros( cs , 4UL );
      checkRows    ( mat, 4UL );
      checkColumns ( mat, 5UL );
      checkNonZeros( mat, 4UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) != 12 || cs(1,1) != 2 ||
          cs(2,0) != 10 || cs(2,1) != 0 ||
          cs(3,0) != -6 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n( 12  2 )\n( 10  0 )\n( -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 || mat(0,3) !=  0 || mat(0,4) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 2 || mat(1,2) != 0 || mat(1,3) != 12 || mat(1,4) != 0 ||
          mat(2,0) != 0 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) != 10 || mat(2,4) != 0 ||
          mat(3,0) != 0 || mat(3,1) != 0 || mat(3,2) != 0 || mat(3,3) != -6 || mat(3,4) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0  0  0  0  0 )\n"
                                     "( 0  2  0 12  0 )\n"
                                     "( 0  0  0 10  0 )\n"
                                     "( 0  0  0 -6  0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Columns Schur product assignment (aliasing)";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 4UL } );
      cs %= blaze::columns( tmat_, { 2UL, 3UL } );

      checkRows    ( cs   , 4UL );
      checkColumns ( cs   , 2UL );
      checkNonZeros( cs   , 4UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 7UL );

      if( cs(0,0) !=   0 || cs(0,1) !=   0 ||
          cs(1,0) !=   0 || cs(1,1) != -32 ||
          cs(2,0) != -15 || cs(2,1) !=  45 ||
          cs(3,0) !=   0 || cs(3,1) != -60 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(   0   0 )\n(   0 -32 )\n( -15  45 )\n(   0 -60 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=   0 || tmat_(0,4) !=   0 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 1 || tmat_(1,2) !=  0 || tmat_(1,3) !=   0 || tmat_(1,4) != -32 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) != -15 || tmat_(2,4) !=  45 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) !=   0 || tmat_(3,4) != -60 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0   0  -2   0   0 )\n"
                                     "( 0   1   0   0  32 )\n"
                                     "( 0   0  -3 -15  45 )\n"
                                     "( 0   0   0   0 -60 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix Schur product assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ {  0, 0 },
                                                           { -1, 2 },
                                                           {  0, 1 },
                                                           { -2, 0 } };

      cs %= mat;

      checkRows    ( cs   , 4UL );
      checkColumns ( cs   , 2UL );
      checkNonZeros( cs   , 3UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) != -4 || cs(1,1) != 2 ||
          cs(2,0) !=  0 || cs(2,1) != 0 ||
          cs(3,0) != 12 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n( -4  2 )\n(  0  0 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 2 || tmat_(1,2) !=  0 || tmat_(1,3) != -4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  2  0 -4 -8 )\n"
                                     "( 0  0 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix Schur product assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ {  0, 0 },
                                                              { -1, 2 },
                                                              {  0, 1 },
                                                              { -2, 0 } };

      cs %= mat;

      checkRows    ( cs   , 4UL );
      checkColumns ( cs   , 2UL );
      checkNonZeros( cs   , 3UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) != -4 || cs(1,1) != 2 ||
          cs(2,0) !=  0 || cs(2,1) != 0 ||
          cs(3,0) != 12 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n( -4  2 )\n(  0  0 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 2 || tmat_(1,2) !=  0 || tmat_(1,3) != -4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  2  0 -4 -8 )\n"
                                     "( 0  0 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix Schur product assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix Schur product assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat{ {  0, 0 },
                                                              { -1, 2 },
                                                              {  0, 1 },
                                                              { -2, 0 } };

      cs %= mat;

      checkRows    ( cs   , 4UL );
      checkColumns ( cs   , 2UL );
      checkNonZeros( cs   , 3UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) != -4 || cs(1,1) != 2 ||
          cs(2,0) !=  0 || cs(2,1) != 0 ||
          cs(3,0) != 12 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n( -4  2 )\n(  0  0 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 2 || tmat_(1,2) !=  0 || tmat_(1,3) != -4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  2  0 -4 -8 )\n"
                                     "( 0  0 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix Schur product assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat{ {  0, 0 },
                                                                 { -1, 2 },
                                                                 {  0, 1 },
                                                                 { -2, 0 } };

      cs %= mat;

      checkRows    ( cs   , 4UL );
      checkColumns ( cs   , 2UL );
      checkNonZeros( cs   , 3UL );
      checkRows    ( tmat_, 4UL );
      checkColumns ( tmat_, 5UL );
      checkNonZeros( tmat_, 9UL );

      if( cs(0,0) !=  0 || cs(0,1) != 0 ||
          cs(1,0) != -4 || cs(1,1) != 2 ||
          cs(2,0) !=  0 || cs(2,1) != 0 ||
          cs(3,0) != 12 || cs(3,1) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(  0  0 )\n( -4  2 )\n(  0  0 )\n( 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != 0 || tmat_(0,1) != 0 || tmat_(0,2) != -2 || tmat_(0,3) !=  0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != 0 || tmat_(1,1) != 2 || tmat_(1,2) !=  0 || tmat_(1,3) != -4 || tmat_(1,4) != -8 ||
          tmat_(2,0) != 0 || tmat_(2,1) != 0 || tmat_(2,2) != -3 || tmat_(2,3) !=  0 || tmat_(2,4) !=  9 ||
          tmat_(3,0) != 0 || tmat_(3,1) != 0 || tmat_(3,2) !=  0 || tmat_(3,3) != 12 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Schur product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( 0  0 -2  0  7 )\n"
                                     "( 0  2  0 -4 -8 )\n"
                                     "( 0  0 -3  0  9 )\n"
                                     "( 0  0  0 12 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Columns multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the Columns
// specialization. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void SparseGeneralTest::testMultAssign()
{
   //=====================================================================================
   // Row-major Columns multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major Columns multiplication assignment (no aliasing)";

      initialize();

      MT mat{ { 0,  0, -2,  0,  7 },
              { 0,  1,  0,  4, -8 },
              { 0,  0, -3,  5,  9 },
              { 0,  0,  0, -6, 10 } };

      auto cs = blaze::columns( mat, { 2UL, 0UL, 3UL, 1UL } );
      cs *= blaze::columns( mat_, { 1UL, 2UL, 2UL, 1UL } );

      checkRows    ( cs ,  4UL );
      checkColumns ( cs ,  4UL );
      checkNonZeros( cs ,  8UL );
      checkRows    ( mat,  4UL );
      checkColumns ( mat,  5UL );
      checkNonZeros( mat, 12UL );

      if( cs(0,0) != 0 || cs(0,1) !=   4 || cs(0,2) !=   4 || cs(0,3) != 0 ||
          cs(1,0) != 0 || cs(1,1) != -12 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != 0 || cs(2,1) !=  -9 || cs(2,2) !=  -9 || cs(2,3) != 0 ||
          cs(3,0) != 0 || cs(3,1) !=  18 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 0   4   4  0 )\n"
                                     "( 0 -12 -12  0 )\n"
                                     "( 0  -9  -9  0 )\n"
                                     "( 0  18  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) !=   4 || mat(0,1) != 0 || mat(0,2) != 0 || mat(0,3) !=   4 || mat(0,4) !=  7 ||
          mat(1,0) != -12 || mat(1,1) != 0 || mat(1,2) != 0 || mat(1,3) != -12 || mat(1,4) != -8 ||
          mat(2,0) !=  -9 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) !=  -9 || mat(2,4) !=  9 ||
          mat(3,0) !=  18 || mat(3,1) != 0 || mat(3,2) != 0 || mat(3,3) !=  18 || mat(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(   4  0  0   4  7 )\n"
                                     "( -12  0  0 -12 -8 )\n"
                                     "(  -9  0  0  -9  9 )\n"
                                     "(  18  0  0  18 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major Columns multiplication assignment (aliasing)";

      initialize();

      auto cs = blaze::columns( mat_, { 2UL, 0UL, 3UL, 1UL } );
      cs *= blaze::columns( mat_, { 1UL, 2UL, 2UL, 1UL } );

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  4UL );
      checkNonZeros( cs  ,  8UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 12UL );

      if( cs(0,0) != 0 || cs(0,1) !=   4 || cs(0,2) !=   4 || cs(0,3) != 0 ||
          cs(1,0) != 0 || cs(1,1) != -12 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != 0 || cs(2,1) !=  -9 || cs(2,2) !=  -9 || cs(2,3) != 0 ||
          cs(3,0) != 0 || cs(3,1) !=  18 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 0   4   4  0 )\n"
                                     "( 0 -12 -12  0 )\n"
                                     "( 0  -9  -9  0 )\n"
                                     "( 0  18  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) !=   4 || mat_(0,1) != 0 || mat_(0,2) != 0 || mat_(0,3) !=   4 || mat_(0,4) !=  7 ||
          mat_(1,0) != -12 || mat_(1,1) != 0 || mat_(1,2) != 0 || mat_(1,3) != -12 || mat_(1,4) != -8 ||
          mat_(2,0) !=  -9 || mat_(2,1) != 0 || mat_(2,2) != 0 || mat_(2,3) !=  -9 || mat_(2,4) !=  9 ||
          mat_(3,0) !=  18 || mat_(3,1) != 0 || mat_(3,2) != 0 || mat_(3,3) !=  18 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n(   4  0  0   4  7 )\n"
                                     "( -12  0  0 -12 -8 )\n"
                                     "(  -9  0  0  -9  9 )\n"
                                     "(  18  0  0  18 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major dense matrix multiplication assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ {  0,  1,  0,  0 },
                                                           { -2,  0, -3,  0 },
                                                           { -2,  0, -3,  0 },
                                                           {  0,  1,  0,  0 } };

      cs *= mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  4UL );
      checkNonZeros( cs  ,  9UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 13UL );

      if( cs(0,0) !=   0 || cs(0,1) != -2 || cs(0,2) !=   0 || cs(0,3) != 0 ||
          cs(1,0) !=  -8 || cs(1,1) !=  1 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != -10 || cs(2,1) != -3 || cs(2,2) != -15 || cs(2,3) != 0 ||
          cs(3,0) !=  12 || cs(3,1) !=  0 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(   0  -2   0  0 )\n"
                                     "(  -8   1 -12  0 )\n"
                                     "( -10  -3 -15  0 )\n"
                                     "(  12   0  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 || mat_(0,4) !=  7 ||
          mat_(1,0) !=  1 || mat_(1,1) != 0 || mat_(1,2) !=  -8 || mat_(1,3) != -12 || mat_(1,4) != -8 ||
          mat_(2,0) != -3 || mat_(2,1) != 0 || mat_(2,2) != -10 || mat_(2,3) != -15 || mat_(2,4) !=  9 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  12 || mat_(3,3) !=  18 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0   0   0  7 )\n"
                                     "(  1  0  -8 -12 -8 )\n"
                                     "( -3  0 -10 -15  9 )\n"
                                     "(  0  0  12  18 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major dense matrix multiplication assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ {  0,  1,  0,  0 },
                                                              { -2,  0, -3,  0 },
                                                              { -2,  0, -3,  0 },
                                                              {  0,  1,  0,  0 } };

      cs *= mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  4UL );
      checkNonZeros( cs  ,  9UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 13UL );

      if( cs(0,0) !=   0 || cs(0,1) != -2 || cs(0,2) !=   0 || cs(0,3) != 0 ||
          cs(1,0) !=  -8 || cs(1,1) !=  1 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != -10 || cs(2,1) != -3 || cs(2,2) != -15 || cs(2,3) != 0 ||
          cs(3,0) !=  12 || cs(3,1) !=  0 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(   0  -2   0  0 )\n"
                                     "(  -8   1 -12  0 )\n"
                                     "( -10  -3 -15  0 )\n"
                                     "(  12   0  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 || mat_(0,4) !=  7 ||
          mat_(1,0) !=  1 || mat_(1,1) != 0 || mat_(1,2) !=  -8 || mat_(1,3) != -12 || mat_(1,4) != -8 ||
          mat_(2,0) != -3 || mat_(2,1) != 0 || mat_(2,2) != -10 || mat_(2,3) != -15 || mat_(2,4) !=  9 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  12 || mat_(3,3) !=  18 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0   0   0  7 )\n"
                                     "(  1  0  -8 -12 -8 )\n"
                                     "( -3  0 -10 -15  9 )\n"
                                     "(  0  0  12  18 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major sparse matrix multiplication assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat{ {  0,  1,  0,  0 },
                                                              { -2,  0, -3,  0 },
                                                              { -2,  0, -3,  0 },
                                                              {  0,  1,  0,  0 } };

      cs *= mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  4UL );
      checkNonZeros( cs  ,  9UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 13UL );

      if( cs(0,0) !=   0 || cs(0,1) != -2 || cs(0,2) !=   0 || cs(0,3) != 0 ||
          cs(1,0) !=  -8 || cs(1,1) !=  1 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != -10 || cs(2,1) != -3 || cs(2,2) != -15 || cs(2,3) != 0 ||
          cs(3,0) !=  12 || cs(3,1) !=  0 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(   0  -2   0  0 )\n"
                                     "(  -8   1 -12  0 )\n"
                                     "( -10  -3 -15  0 )\n"
                                     "(  12   0  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 || mat_(0,4) !=  7 ||
          mat_(1,0) !=  1 || mat_(1,1) != 0 || mat_(1,2) !=  -8 || mat_(1,3) != -12 || mat_(1,4) != -8 ||
          mat_(2,0) != -3 || mat_(2,1) != 0 || mat_(2,2) != -10 || mat_(2,3) != -15 || mat_(2,4) !=  9 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  12 || mat_(3,3) !=  18 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0   0   0  7 )\n"
                                     "(  1  0  -8 -12 -8 )\n"
                                     "( -3  0 -10 -15  9 )\n"
                                     "(  0  0  12  18 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major sparse matrix multiplication assignment";

      initialize();

      auto cs = blaze::columns( mat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat{ {  0,  1,  0,  0 },
                                                                 { -2,  0, -3,  0 },
                                                                 { -2,  0, -3,  0 },
                                                                 {  0,  1,  0,  0 } };

      cs *= mat;

      checkRows    ( cs  ,  4UL );
      checkColumns ( cs  ,  4UL );
      checkNonZeros( cs  ,  9UL );
      checkRows    ( mat_,  4UL );
      checkColumns ( mat_,  5UL );
      checkNonZeros( mat_, 13UL );

      if( cs(0,0) !=   0 || cs(0,1) != -2 || cs(0,2) !=   0 || cs(0,3) != 0 ||
          cs(1,0) !=  -8 || cs(1,1) !=  1 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != -10 || cs(2,1) != -3 || cs(2,2) != -15 || cs(2,3) != 0 ||
          cs(3,0) !=  12 || cs(3,1) !=  0 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(   0  -2   0  0 )\n"
                                     "(  -8   1 -12  0 )\n"
                                     "( -10  -3 -15  0 )\n"
                                     "(  12   0  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat_(0,0) != -2 || mat_(0,1) != 0 || mat_(0,2) !=   0 || mat_(0,3) !=   0 || mat_(0,4) !=  7 ||
          mat_(1,0) !=  1 || mat_(1,1) != 0 || mat_(1,2) !=  -8 || mat_(1,3) != -12 || mat_(1,4) != -8 ||
          mat_(2,0) != -3 || mat_(2,1) != 0 || mat_(2,2) != -10 || mat_(2,3) != -15 || mat_(2,4) !=  9 ||
          mat_(3,0) !=  0 || mat_(3,1) != 0 || mat_(3,2) !=  12 || mat_(3,3) !=  18 || mat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat_ << "\n"
             << "   Expected result:\n( -2  0   0   0  7 )\n"
                                     "(  1  0  -8 -12 -8 )\n"
                                     "( -3  0 -10 -15  9 )\n"
                                     "(  0  0  12  18 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major Columns multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major Columns multiplication assignment (no aliasing)";

      initialize();

      OMT mat{ { 0,  0, -2,  0,  7 },
               { 0,  1,  0,  4, -8 },
               { 0,  0, -3,  5,  9 },
               { 0,  0,  0, -6, 10 } };

      auto cs = blaze::columns( mat, { 2UL, 0UL, 3UL, 1UL } );
      cs *= blaze::columns( tmat_, { 1UL, 2UL, 2UL, 1UL } );

      checkRows    ( cs ,  4UL );
      checkColumns ( cs ,  4UL );
      checkNonZeros( cs ,  8UL );
      checkRows    ( mat,  4UL );
      checkColumns ( mat,  5UL );
      checkNonZeros( mat, 12UL );

      if( cs(0,0) != 0 || cs(0,1) !=   4 || cs(0,2) !=   4 || cs(0,3) != 0 ||
          cs(1,0) != 0 || cs(1,1) != -12 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != 0 || cs(2,1) !=  -9 || cs(2,2) !=  -9 || cs(2,3) != 0 ||
          cs(3,0) != 0 || cs(3,1) !=  18 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 0   4   4  0 )\n"
                                     "( 0 -12 -12  0 )\n"
                                     "( 0  -9  -9  0 )\n"
                                     "( 0  18  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( mat(0,0) !=   4 || mat(0,1) != 0 || mat(0,2) != 0 || mat(0,3) !=   4 || mat(0,4) !=  7 ||
          mat(1,0) != -12 || mat(1,1) != 0 || mat(1,2) != 0 || mat(1,3) != -12 || mat(1,4) != -8 ||
          mat(2,0) !=  -9 || mat(2,1) != 0 || mat(2,2) != 0 || mat(2,3) !=  -9 || mat(2,4) !=  9 ||
          mat(3,0) !=  18 || mat(3,1) != 0 || mat(3,2) != 0 || mat(3,3) !=  18 || mat(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n(   4  0  0   4  7 )\n"
                                     "( -12  0  0 -12 -8 )\n"
                                     "(  -9  0  0  -9  9 )\n"
                                     "(  18  0  0  18 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major Columns multiplication assignment (aliasing)";

      initialize();

      auto cs = blaze::columns( tmat_, { 2UL, 0UL, 3UL, 1UL } );
      cs *= blaze::columns( tmat_, { 1UL, 2UL, 2UL, 1UL } );

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  4UL );
      checkNonZeros( cs   ,  8UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 12UL );

      if( cs(0,0) != 0 || cs(0,1) !=   4 || cs(0,2) !=   4 || cs(0,3) != 0 ||
          cs(1,0) != 0 || cs(1,1) != -12 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != 0 || cs(2,1) !=  -9 || cs(2,2) !=  -9 || cs(2,3) != 0 ||
          cs(3,0) != 0 || cs(3,1) !=  18 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n( 0   4   4  0 )\n"
                                     "( 0 -12 -12  0 )\n"
                                     "( 0  -9  -9  0 )\n"
                                     "( 0  18  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) !=   4 || tmat_(0,1) != 0 || tmat_(0,2) != 0 || tmat_(0,3) !=   4 || tmat_(0,4) !=  7 ||
          tmat_(1,0) != -12 || tmat_(1,1) != 0 || tmat_(1,2) != 0 || tmat_(1,3) != -12 || tmat_(1,4) != -8 ||
          tmat_(2,0) !=  -9 || tmat_(2,1) != 0 || tmat_(2,2) != 0 || tmat_(2,3) !=  -9 || tmat_(2,4) !=  9 ||
          tmat_(3,0) !=  18 || tmat_(3,1) != 0 || tmat_(3,2) != 0 || tmat_(3,3) !=  18 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n(   4  0  0   4  7 )\n"
                                     "( -12  0  0 -12 -8 )\n"
                                     "(  -9  0  0  -9  9 )\n"
                                     "(  18  0  0  18 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major dense matrix multiplication assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat{ {  0,  1,  0,  0 },
                                                           { -2,  0, -3,  0 },
                                                           { -2,  0, -3,  0 },
                                                           {  0,  1,  0,  0 } };

      cs *= mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  4UL );
      checkNonZeros( cs   ,  9UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 13UL );

      if( cs(0,0) !=   0 || cs(0,1) != -2 || cs(0,2) !=   0 || cs(0,3) != 0 ||
          cs(1,0) !=  -8 || cs(1,1) !=  1 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != -10 || cs(2,1) != -3 || cs(2,2) != -15 || cs(2,3) != 0 ||
          cs(3,0) !=  12 || cs(3,1) !=  0 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(   0  -2   0  0 )\n"
                                     "(  -8   1 -12  0 )\n"
                                     "( -10  -3 -15  0 )\n"
                                     "(  12   0  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) != 0 || tmat_(0,2) !=   0 || tmat_(0,3) !=   0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) !=  1 || tmat_(1,1) != 0 || tmat_(1,2) !=  -8 || tmat_(1,3) != -12 || tmat_(1,4) != -8 ||
          tmat_(2,0) != -3 || tmat_(2,1) != 0 || tmat_(2,2) != -10 || tmat_(2,3) != -15 || tmat_(2,4) !=  9 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != 0 || tmat_(3,2) !=  12 || tmat_(3,3) !=  18 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0   0   0  7 )\n"
                                     "(  1  0  -8 -12 -8 )\n"
                                     "( -3  0 -10 -15  9 )\n"
                                     "(  0  0  12  18 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major dense matrix multiplication assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat{ {  0,  1,  0,  0 },
                                                              { -2,  0, -3,  0 },
                                                              { -2,  0, -3,  0 },
                                                              {  0,  1,  0,  0 } };

      cs *= mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  4UL );
      checkNonZeros( cs   ,  9UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 13UL );

      if( cs(0,0) !=   0 || cs(0,1) != -2 || cs(0,2) !=   0 || cs(0,3) != 0 ||
          cs(1,0) !=  -8 || cs(1,1) !=  1 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != -10 || cs(2,1) != -3 || cs(2,2) != -15 || cs(2,3) != 0 ||
          cs(3,0) !=  12 || cs(3,1) !=  0 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(   0  -2   0  0 )\n"
                                     "(  -8   1 -12  0 )\n"
                                     "( -10  -3 -15  0 )\n"
                                     "(  12   0  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) != 0 || tmat_(0,2) !=   0 || tmat_(0,3) !=   0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) !=  1 || tmat_(1,1) != 0 || tmat_(1,2) !=  -8 || tmat_(1,3) != -12 || tmat_(1,4) != -8 ||
          tmat_(2,0) != -3 || tmat_(2,1) != 0 || tmat_(2,2) != -10 || tmat_(2,3) != -15 || tmat_(2,4) !=  9 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != 0 || tmat_(3,2) !=  12 || tmat_(3,3) !=  18 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0   0   0  7 )\n"
                                     "(  1  0  -8 -12 -8 )\n"
                                     "( -3  0 -10 -15  9 )\n"
                                     "(  0  0  12  18 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix multiplication assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major sparse matrix multiplication assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::rowMajor> mat{ {  0,  1,  0,  0 },
                                                              { -2,  0, -3,  0 },
                                                              { -2,  0, -3,  0 },
                                                              {  0,  1,  0,  0 } };

      cs *= mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  4UL );
      checkNonZeros( cs   ,  9UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 13UL );

      if( cs(0,0) !=   0 || cs(0,1) != -2 || cs(0,2) !=   0 || cs(0,3) != 0 ||
          cs(1,0) !=  -8 || cs(1,1) !=  1 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != -10 || cs(2,1) != -3 || cs(2,2) != -15 || cs(2,3) != 0 ||
          cs(3,0) !=  12 || cs(3,1) !=  0 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(   0  -2   0  0 )\n"
                                     "(  -8   1 -12  0 )\n"
                                     "( -10  -3 -15  0 )\n"
                                     "(  12   0  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) != 0 || tmat_(0,2) !=   0 || tmat_(0,3) !=   0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) !=  1 || tmat_(1,1) != 0 || tmat_(1,2) !=  -8 || tmat_(1,3) != -12 || tmat_(1,4) != -8 ||
          tmat_(2,0) != -3 || tmat_(2,1) != 0 || tmat_(2,2) != -10 || tmat_(2,3) != -15 || tmat_(2,4) !=  9 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != 0 || tmat_(3,2) !=  12 || tmat_(3,3) !=  18 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0   0   0  7 )\n"
                                     "(  1  0  -8 -12 -8 )\n"
                                     "( -3  0 -10 -15  9 )\n"
                                     "(  0  0  12  18 10 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major sparse matrix multiplication assignment";

      initialize();

      auto cs = blaze::columns( tmat_, { 2UL, 0UL, 3UL, 1UL } );

      const blaze::CompressedMatrix<int,blaze::columnMajor> mat{ {  0,  1,  0,  0 },
                                                                 { -2,  0, -3,  0 },
                                                                 { -2,  0, -3,  0 },
                                                                 {  0,  1,  0,  0 } };

      cs *= mat;

      checkRows    ( cs   ,  4UL );
      checkColumns ( cs   ,  4UL );
      checkNonZeros( cs   ,  9UL );
      checkRows    ( tmat_,  4UL );
      checkColumns ( tmat_,  5UL );
      checkNonZeros( tmat_, 13UL );

      if( cs(0,0) !=   0 || cs(0,1) != -2 || cs(0,2) !=   0 || cs(0,3) != 0 ||
          cs(1,0) !=  -8 || cs(1,1) !=  1 || cs(1,2) != -12 || cs(1,3) != 0 ||
          cs(2,0) != -10 || cs(2,1) != -3 || cs(2,2) != -15 || cs(2,3) != 0 ||
          cs(3,0) !=  12 || cs(3,1) !=  0 || cs(3,2) !=  18 || cs(3,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << cs << "\n"
             << "   Expected result:\n(   0  -2   0  0 )\n"
                                     "(  -8   1 -12  0 )\n"
                                     "( -10  -3 -15  0 )\n"
                                     "(  12   0  18  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( tmat_(0,0) != -2 || tmat_(0,1) != 0 || tmat_(0,2) !=   0 || tmat_(0,3) !=   0 || tmat_(0,4) !=  7 ||
          tmat_(1,0) !=  1 || tmat_(1,1) != 0 || tmat_(1,2) !=  -8 || tmat_(1,3) != -12 || tmat_(1,4) != -8 ||
          tmat_(2,0) != -3 || tmat_(2,1) != 0 || tmat_(2,2) != -10 || tmat_(2,3) != -15 || tmat_(2,4) !=  9 ||
          tmat_(3,0) !=  0 || tmat_(3,1) != 0 || tmat_(3,2) !=  12 || tmat_(3,3) !=  18 || tmat_(3,4) != 10 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << tmat_ << "\n"
             << "   Expected result:\n( -2  0   0   0  7 )\n"
                                     "(  1  0  -8 -12 -8 )\n"
                                     "( -3  0 -10 -15  9 )\n"
                                     "(  0  0  12  18 10 )\n";
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
void SparseGeneralTest::initialize()
{
   // Initializing the row-major compressed matrix
   mat_.reset();
   mat_(1,1) =  1;
   mat_(0,2) = -2;
   mat_(2,2) = -3;
   mat_(1,3) =  4;
   mat_(2,3) =  5;
   mat_(3,3) = -6;
   mat_(0,4) =  7;
   mat_(1,4) = -8;
   mat_(2,4) =  9;
   mat_(3,4) = 10;

   // Initializing the column-major compressed matrix
   tmat_.reset();
   tmat_(1,1) =  1;
   tmat_(0,2) = -2;
   tmat_(2,2) = -3;
   tmat_(1,3) =  4;
   tmat_(2,3) =  5;
   tmat_(3,3) = -6;
   tmat_(0,4) =  7;
   tmat_(1,4) = -8;
   tmat_(2,4) =  9;
   tmat_(3,4) = 10;
}
//*************************************************************************************************

} // namespace columns

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
   std::cout << "   Running Columns sparse general test (part 1)..." << std::endl;

   try
   {
      RUN_COLUMNS_SPARSEGENERAL_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during Columns sparse general test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
