//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/uniuppermatrix/RowTest.h
//  \brief Header file for the UniUpperMatrix row test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_UNIUPPERMATRIX_ROWTEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_UNIUPPERMATRIX_ROWTEST_H_


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
#include <blaze/math/UniUpperMatrix.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace uniuppermatrix {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for assignment tests to a single row of a UniUpperMatrix.
//
// This class performs assignment tests to a single row of a UniUpperMatrix. It performs a series
// of both compile time as well as runtime tests.
*/
class RowTest
{
 private:
   //**Type definitions****************************************************************************
   //! Type of the dense upper unitriangular matrix.
   using DUT = blaze::UniUpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> >;

   //! Opposite dense upper unitriangular matrix type.
   using DOUT = DUT::OppositeType;

   //! Type of the sparse upper unitriangular matrix.
   using SUT = blaze::UniUpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> >;

   //! Opposite sparse upper unitriangular matrix type.
   using SOUT = SUT::OppositeType;
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
   template< typename UT > void testAssignment();
   template< typename UT > void testAddAssign ();
   template< typename UT > void testSubAssign ();
   template< typename UT > void testMultAssign();

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
   template< typename UT > void init( UT& upper );
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
/*!\brief Test of the assignment to rows of a UniUpperMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to a single row of a UniUpperMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename UT >  // Type of the uniupper matrix
void RowTest::testAssignment()
{
   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  1 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Dense vector assignment test 1";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[1] =  1;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );
      row1 = vec;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 6UL );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0  1 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -2 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  1 -2 )\n( 0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  0 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Dense vector assignment test 2";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 9  1 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Dense vector assignment test 3";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL );
      vec[0] =  9;
      vec[0] =  1;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  1 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Sparse vector assignment test 1";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[1] =  1;
      vec[2] = -2;
      vec.insert( 0UL, 0 );

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );
      row1 = vec;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 6UL );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0  1 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -2 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  1 -2 )\n( 0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  0 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Sparse vector assignment test 2";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 1UL );
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 9  1 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Sparse vector assignment test 3";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[0] =  9;
      vec[0] =  1;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the addition assignment to rows of a UniUpperMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment to a single row of a UniUpperMatrix. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename UT >  // Type of the uniupper matrix
void RowTest::testAddAssign()
{
   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  1 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Dense vector addition assignment test 1";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );
      row1 += vec;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 6UL );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0  1 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -2 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  1 -2 )\n( 0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  2 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Dense vector addition assignment test 2";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[1] =  1;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 9  1 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Dense vector addition assignment test 3";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] =  9;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  1 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Sparse vector addition assignment test 1";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[2] = -2;
      vec.insert( 0UL, 0 );

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );
      row1 += vec;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 6UL );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0  1 -2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -2 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  1 -2 )\n( 0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  2 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Sparse vector addition assignment test 2";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[1] =  1;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 9  1 -2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Sparse vector addition assignment test 3";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] =  9;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the subtraction assignment to rows of a UniUpperMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment to a single row of a UniUpperMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename UT >  // Type of the uniupper matrix
void RowTest::testSubAssign()
{
   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  1  2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Dense vector subtraction assignment test 1";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );
      row1 -= vec;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 6UL );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0  1 2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 2 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  1  2 )\n( 0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  0  2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Dense vector subtraction assignment test 2";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[1] =  1;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 )      (  1 -4  7 )
   // ( 0  1  0 )  =>  ( -9  1  2 )
   // ( 0  0  1 )      (  0  0  1 )
   {
      test_ = "Dense vector subtraction assignment test 3";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );
      vec[0] =  9;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  1  2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Sparse vector subtraction assignment test 1";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[2] = -2;
      vec.insert( 0UL, 0 );

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );
      row1 -= vec;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 6UL );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0  1  2 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 2 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  1  2 )\n( 0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  0  2 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Sparse vector subtraction assignment test 2";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[1] =  1;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 )      (  1 -4  7 )
   // ( 0  1  0 )  =>  ( -9  1  2 )
   // ( 0  0  1 )      (  0  0  1 )
   {
      test_ = "Sparse vector subtraction assignment test 3";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 2UL );
      vec[0] =  9;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication assignment to rows of a UniUpperMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment to a single row of a
// UniUpperMatrix. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename UT >  // Type of the uniupper matrix
void RowTest::testMultAssign()
{
   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  1  0 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Dense vector multiplication assignment test 1";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL );
      vec[0] =  9;
      vec[1] =  1;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );
      row1 *= vec;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 5UL );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0  1  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  1  0 )\n( 0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  0  0 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Dense vector multiplication assignment test 2";

      blaze::DynamicVector<int,blaze::rowVector> vec( 3UL, 0 );

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 *= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  1  0 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Sparse vector multiplication assignment test 1";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL, 3UL );
      vec[0] =  9;
      vec[1] =  1;
      vec[2] = -2;

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );
      row1 *= vec;

      checkRows    ( upper, 3UL );
      checkColumns ( upper, 3UL );
      checkNonZeros( upper, 5UL );

      if( row1[0] != 0 || row1[1] != 1 || row1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << row1 << "\n"
             << "   Expected result:\n( 0  1  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to row failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 )\n( 0  1  0 )\n( 0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 )      ( 1 -4  7 )
   // ( 0  1  0 )  =>  ( 0  0  0 )
   // ( 0  0  1 )      ( 0  0  1 )
   {
      test_ = "Sparse vector multiplication assignment test 2";

      blaze::CompressedVector<int,blaze::rowVector> vec( 3UL );

      UT upper;
      init( upper );

      auto row1 = row( upper, 1UL );

      try {
         row1 *= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
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
/*!\brief Initializing the given upper unitriangular matrix.
//
// \return void
//
// This function is called before each test case to initialize the given upper unitriangular
// matrix.
*/
template< typename UT >
void RowTest::init( UT& upper )
{
   upper.resize( 3UL );
   upper(0,1) = -4;
   upper(0,2) =  7;
   upper(1,2) =  0;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the assignment to a single row of a UniUpperMatrix.
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
/*!\brief Macro for the execution of the UniUpperMatrix row test.
*/
#define RUN_UNIUPPERMATRIX_ROW_TEST \
   blazetest::mathtest::adaptors::uniuppermatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace uniuppermatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
