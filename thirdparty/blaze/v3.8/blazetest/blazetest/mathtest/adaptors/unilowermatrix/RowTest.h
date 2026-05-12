//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/unilowermatrix/RowTest.h
//  \brief Header file for the UniLowerMatrix row test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_UNILOWERMATRIX_ROWTEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_UNILOWERMATRIX_ROWTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Row.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace unilowermatrix {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for assignment tests to a single row of a UniLowerMatrix.
//
// This class performs assignment tests to a single row of a UniLowerMatrix. It performs a series
// of both compile time as well as runtime tests.
*/
class RowTest
{
 private:
   //**Type definitions****************************************************************************
   //! Type of the dense lower unitriangular matrix.
   using DLT = blaze::UniLowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> >;

   //! Opposite dense lower unitriangular matrix type.
   using DOLT = DLT::OppositeType;

   //! Type of the sparse lower unitriangular matrix.
   using SLT = blaze::UniLowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> >;

   //! Opposite sparse lower unitriangular matrix type.
   using SOLT = SLT::OppositeType;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit RowTest();
   // No explicitly declared copy constructor.
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

 private:
   //**Test functions******************************************************************************
   /*!\name Test functions */
   //@{
   template< typename LT > void testAssignment();
   template< typename LT > void testAddAssign ();
   template< typename LT > void testSubAssign ();
   template< typename LT > void testMultAssign();

   template< typename Type >
   void checkRows( const Type& matrix, size_t expectedRows ) const;

   template< typename Type >
   void checkColumns( const Type& matrix, size_t expectedColumns ) const;

   template< typename Type >
   void checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< typename LT > void init( LT& lower );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the assignment to rows of a UniLowerMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to a single row of a UniLowerMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename LT >  // Type of the unilower matrix
void RowTest::testAssignment()
{
   //=====================================================================================
   // Row-major dense vector assignment
   //=====================================================================================

   // (  1  0  0 )      ( 1  0  0 )
   // ( -4  1  0 )  =>  ( 2  1  0 )
   // (  7  0  1 )      ( 7  0  1 )
   {
      test_ = "Dense vector assignment test 1";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] = 2;
      vec[1] = 1;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );
      row1 = vec;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );

      if( row1[0] != 2 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 2 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0 )      ( 1  0  0 )
   // ( -4  1  0 )  =>  ( 2  0  0 )
   // (  7  0  1 )      ( 7  0  1 )
   {
      test_ = "Dense vector assignment test 2";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] = 2;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0 )      ( 1  0  0 )
   // ( -4  1  0 )  =>  ( 2  1  9 )
   // (  7  0  1 )      ( 7  0  1 )
   {
      test_ = "Dense vector assignment test 3";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL );
      vec[0] = 2;
      vec[1] = 1;
      vec[2] = 9;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse vector assignment
   //=====================================================================================

   // (  1  0  0 )      ( 1  0  0 )
   // ( -4  1  0 )  =>  ( 2  1  0 )
   // (  7  0  1 )      ( 7  0  1 )
   {
      test_ = "Sparse vector assignment test 1";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[0] = 2;
      vec[1] = 1;
      vec.insert( 2UL, 0 );

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );
      row1 = vec;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );

      if( row1[0] != 2 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 2 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) != 1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != 2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) != 7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 2 1 0 )\n( 7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0 )      ( 1  0  0 )
   // ( -4  1  0 )  =>  ( 2  0  0 )
   // (  7  0  1 )      ( 7  0  1 )
   {
      test_ = "Sparse vector assignment test 2";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 1UL );
      vec[0] = 2;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0 )      ( 1  0  0 )
   // ( -4  1  0 )  =>  ( 2  1  9 )
   // (  7  0  1 )      ( 7  0  1 )
   {
      test_ = "Sparse vector assignment test 3";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[0] = 2;
      vec[1] = 1;
      vec[2] = 9;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the addition assignment to rows of a UniLowerMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment to a single row of a UniLowerMatrix. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename LT >  // Type of the unilower matrix
void RowTest::testAddAssign()
{
   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -2  1  0 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Dense vector addition assignment test 1";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] = 2;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );
      row1 += vec;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );

      if( row1[0] != -2 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -2 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -2  2  0 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Dense vector addition assignment test 2";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] = 2;
      vec[1] = 1;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -2  1  9 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Dense vector addition assignment test 3";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] = 2;
      vec[2] = 9;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -2  1  0 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Sparse vector addition assignment test 1";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] = 2;
      vec.insert( 1UL, 0 );

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );
      row1 += vec;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );

      if( row1[0] != -2 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -2 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -2 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -2 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -2  2  0 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Sparse vector addition assignment test 2";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] = 2;
      vec[1] = 1;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -2  1  9 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Sparse vector addition assignment test 3";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] = 2;
      vec[2] = 9;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the subtraction assignment to rows of a UniLowerMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment to a single row of a UniLowerMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename LT >  // Type of the unilower matrix
void RowTest::testSubAssign()
{
   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -6  1  0 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Dense vector subtraction assignment test 1";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] = 2;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );
      row1 -= vec;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );

      if( row1[0] != -6 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -6 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -6 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -6 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -6  0  0 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Dense vector subtraction assignment test 2";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] = 2;
      vec[1] = 1;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -6  1  9 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Dense vector subtraction assignment test 3";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] = 2;
      vec[2] = 9;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -6  1  0 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Sparse vector subtraction assignment test 1";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] = 2;
      vec.insert( 1UL, 0 );

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );
      row1 -= vec;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );

      if( row1[0] != -6 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -6 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -6 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -6 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -6  0  0 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Sparse vector subtraction assignment test 2";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] = 2;
      vec[1] = 1;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -6  1  9 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Sparse vector subtraction assignment test 3";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] = 2;
      vec[2] = 9;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication assignment to rows of a UniLowerMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment to a single row of a
// UniLowerMatrix. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename LT >  // Type of the unilower matrix
void RowTest::testMultAssign()
{
   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -8  1  0 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Dense vector multiplication assignment test 1";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL );
      vec[0] = 2;
      vec[1] = 1;
      vec[2] = 3;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );
      row1 *= vec;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );

      if( row1[0] != -8 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -8 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -8 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0 )      ( 1  0  0 )
   // ( -4  1  0 )  =>  ( 0  0  0 )
   // (  7  0  1 )      ( 7  0  1 )
   {
      test_ = "Dense vector multiplication assignment test 2";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 *= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   // (  1  0  0 )      (  1  0  0 )
   // ( -4  1  0 )  =>  ( -8  1  0 )
   // (  7  0  1 )      (  7  0  1 )
   {
      test_ = "Sparse vector multiplication assignment test 1";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[0] = 2;
      vec[1] = 1;
      vec[2] = 3;

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );
      row1 *= vec;

      checkRows    ( lower, 3UL );
      checkColumns ( lower, 3UL );
      checkNonZeros( lower, 5UL );

      if( row1[0] != -8 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( -8 1 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) != 0 || lower(0,2) != 0 ||
          lower(1,0) != -8 || lower(1,1) != 1 || lower(1,2) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 0 || lower(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1 0 0 )\n( -8 1 0 )\n(  7 0 1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0 )      ( 1  0  0 )
   // ( -4  1  0 )  =>  ( 0  0  0 )
   // (  7  0  1 )      ( 7  0  1 )
   {
      test_ = "Sparse vector multiplication assignment test 2";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL );

      LT lower;
      init( lower );

      auto row1 = row( lower, 1UL );

      try {
         row1 *= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of rows of the given matrix.
//
// \param matrix The matrix to be checked.
// \param expectedRows The expected number of rows of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of rows of the given matrix. In case the actual number of
// rows does not correspond to the given expected number of rows, a \a std::runtime_error
// exception is thrown.
*/
template< typename Type >  // Type of the matrix
void RowTest::checkRows( const Type& matrix, size_t expectedRows ) const
{
   if( matrix.rows() != expectedRows ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of rows detected\n"
          << " Details:\n"
          << "   Number of rows         : " << matrix.rows() << "\n"
          << "   Expected number of rows: " << expectedRows << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of columns of the given matrix.
//
// \param matrix The matrix to be checked.
// \param expectedColumns The expected number of columns of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of columns of the given matrix. In case the actual number of
// columns does not correspond to the given expected number of columns, a \a std::runtime_error
// exception is thrown.
*/
template< typename Type >  // Type of the matrix
void RowTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
{
   if( matrix.columns() != expectedColumns ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of columns detected\n"
          << " Details:\n"
          << "   Number of columns         : " << matrix.columns() << "\n"
          << "   Expected number of columns: " << expectedColumns << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of non-zero elements of the given matrix.
//
// \param matrix The matrix to be checked.
// \param expectedNonZeros The expected number of non-zero elements of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given matrix. In case the
// actual number of non-zero elements does not correspond to the given expected number,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the matrix
void RowTest::checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const
{
   if( nonZeros( matrix ) != expectedNonZeros ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of non-zero elements\n"
          << " Details:\n"
          << "   Number of non-zeros         : " << nonZeros( matrix ) << "\n"
          << "   Expected number of non-zeros: " << expectedNonZeros << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( capacity( matrix ) < nonZeros( matrix ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Number of non-zeros: " << nonZeros( matrix ) << "\n"
          << "   Capacity           : " << capacity( matrix ) << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initializing the given lower unitriangular matrix.
//
// \return void
//
// This function is called before each test case to initialize the given lower unitriangular
// matrix.
*/
template< typename LT >
void RowTest::init( LT& lower )
{
   lower.resize( 3UL );
   lower(1,0) = -4;
   lower(2,0) =  7;
   lower(2,1) =  0;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the assignment to a single row of a UniLowerMatrix.
//
// \return void
*/
void runTest()
{
   RowTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the UniLowerMatrix row test.
*/
#define RUN_UNILOWERMATRIX_ROW_TEST \
   blazetest::mathtest::adaptors::unilowermatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace unilowermatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
