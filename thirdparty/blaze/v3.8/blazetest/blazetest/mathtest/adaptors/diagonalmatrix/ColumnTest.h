//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/diagonalmatrix/ColumnTest.h
//  \brief Header file for the DiagonalMatrix column test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_DIAGONALMATRIX_COLUMNTEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_DIAGONALMATRIX_COLUMNTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/Column.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace diagonalmatrix {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for assignment tests to a single column of a DiagonalMatrix.
//
// This class performs assignment tests to a single column of a DiagonalMatrix. It performs a series
// of both compile time as well as runtime tests.
*/
class ColumnTest
{
 private:
   //**Type definitions****************************************************************************
   //! Type of the dense diagonal matrix.
   using DDT = blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> >;

   //! Opposite dense diagonal matrix type.
   using DODT = DDT::OppositeType;

   //! Type of the sparse diagonal matrix.
   using SDT = blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> >;

   //! Opposite sparse diagonal matrix type.
   using SODT = SDT::OppositeType;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit ColumnTest();
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
   template< typename DT > void testAssignment();
   template< typename DT > void testAddAssign ();
   template< typename DT > void testSubAssign ();
   template< typename DT > void testMultAssign();

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
   template< typename DT > void init( DT& diag );
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
/*!\brief Test of the assignment to columns of a DiagonalMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to a single column of a DiagonalMatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename DT >  // Type of the diagonal matrix
void ColumnTest::testAssignment()
{
   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Dense vector assignment test 1";

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[1] = 8;

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );
      col1 = vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  9  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Dense vector assignment test 2";

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[0] = 9;
      vec[1] = 8;

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  2  3 )
   {
      test_ = "Dense vector assignment test 3";

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[1] = 8;
      vec[2] = 2;

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Sparse vector assignment test 1";

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[1] = 8;
      vec.insert( 0UL, 0 );
      vec.insert( 2UL, 0 );

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );
      col1 = vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  9  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Sparse vector assignment test 2";

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = 9;
      vec[1] = 8;
      vec.insert( 2UL, 0 );

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  2  3 )
   {
      test_ = "Sparse vector assignment test 3";

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[1] = 8;
      vec[2] = 2;
      vec.insert( 0UL, 0 );

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the addition assignment to columns of a DiagonalMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment to a single column of a DiagonalMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename DT >  // Type of the diagonal matrix
void ColumnTest::testAddAssign()
{
   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Dense vector addition assignment test 1";

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[1] = 6;

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );
      col1 += vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  9  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Dense vector addition assignment test 2";

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[0] = 9;
      vec[1] = 6;

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  2  3 )
   {
      test_ = "Dense vector addition assignment test 3";

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[1] = 6;
      vec[2] = 2;

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Sparse vector addition assignment test 1";

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[1] = 6;
      vec.insert( 0UL, 0 );
      vec.insert( 2UL, 0 );

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );
      col1 += vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  9  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Sparse vector addition assignment test 2";

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = 9;
      vec[1] = 6;
      vec.insert( 2UL, 0 );

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  2  3 )
   {
      test_ = "Sparse vector addition assignment test 3";

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[1] = 6;
      vec[2] = 2;
      vec.insert( 0UL, 0 );

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the subtraction assignment to columns of a DiagonalMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment to a single column of a
// DiagonalMatrix. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename DT >  // Type of the diagonal matrix
void ColumnTest::testSubAssign()
{
   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Dense vector subtraction assignment test 1";

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[1] = -6;

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );
      col1 -= vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  9  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Dense vector subtraction assignment test 2";

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[0] = -9;
      vec[1] = -6;

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  2  3 )
   {
      test_ = "Dense vector subtraction assignment test 3";

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL, 0 );
      vec[1] = -6;
      vec[2] = -2;

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Sparse vector subtraction assignment test 1";

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[1] = -6;
      vec.insert( 0UL, 0 );
      vec.insert( 2UL, 0 );

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );
      col1 -= vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0 )      ( 1  9  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Sparse vector subtraction assignment test 2";

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = -9;
      vec[1] = -6;
      vec.insert( 2UL, 0 );

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  2  3 )
   {
      test_ = "Sparse vector subtraction assignment test 3";

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[1] = -6;
      vec[2] = -2;
      vec.insert( 0UL, 0 );

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );

      try {
         col1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication assignment to columns of a DiagonalMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment to a single column of a
// DiagonalMatrix. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename DT >  // Type of the diagonal matrix
void ColumnTest::testMultAssign()
{
   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Dense vector multiplication assignment test";

      blaze::DynamicVector<int,blaze::columnVector> vec( 3UL );
      vec[0] = 9;
      vec[1] = 4;
      vec[2] = 2;

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );
      col1 *= vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   // ( 1  0  0 )      ( 1  0  0 )
   // ( 0  2  0 )  =>  ( 0  8  0 )
   // ( 0  0  3 )      ( 0  0  3 )
   {
      test_ = "Sparse vector multiplication assignment test";

      blaze::CompressedVector<int,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = 9;
      vec[1] = 4;
      vec[2] = 2;

      DT diag;
      init( diag );

      auto col1 = column( diag, 1UL );
      col1 *= vec;

      checkRows    ( diag, 3UL );
      checkColumns ( diag, 3UL );
      checkNonZeros( diag, 3UL );

      if( col1[0] != 0 || col1[1] != 8 || col1[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( 0 8 0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) != 0 || diag(0,2) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 8 || diag(1,2) != 0 ||
          diag(2,0) != 0 || diag(2,1) != 0 || diag(2,2) != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 0 8 0 )\n( 0 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
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
void ColumnTest::checkRows( const Type& matrix, size_t expectedRows ) const
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
void ColumnTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
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
void ColumnTest::checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const
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
/*!\brief Initializing the given diagonal matrix.
//
// \return void
//
// This function is called before each test case to initialize the given diagonal matrix.
*/
template< typename DT >
void ColumnTest::init( DT& diag )
{
   diag.resize( 3UL );
   diag(0,0) = 1;
   diag(1,1) = 2;
   diag(2,2) = 3;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the assignment to a single column of a DiagonalMatrix.
//
// \return void
*/
void runTest()
{
   ColumnTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the DiagonalMatrix column test.
*/
#define RUN_DIAGONALMATRIX_COLUMN_TEST \
   blazetest::mathtest::adaptors::diagonalmatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace diagonalmatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
