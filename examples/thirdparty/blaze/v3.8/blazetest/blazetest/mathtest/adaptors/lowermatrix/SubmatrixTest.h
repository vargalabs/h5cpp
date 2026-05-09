//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/lowermatrix/SubmatrixTest.h
//  \brief Header file for the LowerMatrix submatrix test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_LOWERMATRIX_SUBMATRIXTEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_LOWERMATRIX_SUBMATRIXTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/Submatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace lowermatrix {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for assignment tests to a submatrix of a LowerMatrix.
//
// This class performs assignment tests to a submatrix of a LowerMatrix. It performs a series
// of both compile time as well as runtime tests.
*/
class SubmatrixTest
{
 private:
   //**Type definitions****************************************************************************
   //! Type of the dense lower triangular matrix.
   using DLT = blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> >;

   //! Opposite dense lower triangular matrix type.
   using DOLT = DLT::OppositeType;

   //! Type of the sparse lower triangular matrix.
   using SLT = blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> >;

   //! Opposite sparse lower triangular matrix type.
   using SOLT = SLT::OppositeType;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit SubmatrixTest();
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
   template< typename LT > void testAssignment ();
   template< typename LT > void testAddAssign  ();
   template< typename LT > void testSubAssign  ();
   template< typename LT > void testSchurAssign();

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
/*!\brief Test of the assignment to a submatrix of a LowerMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to a submatrix of a LowerMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename LT >  // Type of the lower matrix
void SubmatrixTest::testAssignment()
{
   //=====================================================================================
   // Dense matrix assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 18  0  0 )
   // (  7  0  3  0 )      (  7 14 11  0 )
   // ( -2  0  1  5 )      ( -2 15 19  5 )
   {
      test_ = "Row-major dense matrix assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) != 14 || sm(2,1) != 11 ||
          sm(3,0) != 15 || sm(3,1) != 19 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n( 14 11 )\n( 15 19 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 18 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 14 || lower(2,2) != 11 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 15 || lower(3,2) != 19 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 18  0  0 )\n"
                                     "(  7 14 11  0 )\n"
                                     "( -2 15 19  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( 12 18  0  0 )
   // (  7  0  3  0 )      ( 15 17 11  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major dense matrix assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 15 || sm(1,1) != 17 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 18  0  0 )\n( 15 17 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != 12 || lower(1,1) != 18 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 15 || lower(2,1) != 17 || lower(2,2) != 11 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( 12 18  0  0 )\n"
                                     "( 15 17 11  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 11 12  0 )
   // (  7  0  3  0 )      (  7 13 14  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major dense matrix assignment test 3";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 18  0  0 )
   // (  7  0  3  0 )      (  7 14 11  0 )
   // ( -2  0  1  5 )      ( -2 15 19  5 )
   {
      test_ = "Column-major dense matrix assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) != 14 || sm(2,1) != 11 ||
          sm(3,0) != 15 || sm(3,1) != 19 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n( 14 11 )\n( 15 19 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 18 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 14 || lower(2,2) != 11 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 15 || lower(3,2) != 19 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 18  0  0 )\n"
                                     "(  7 14 11  0 )\n"
                                     "( -2 15 19  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( 12 18  0  0 )
   // (  7  0  3  0 )      ( 15 17 11  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major dense matrix assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 15 || sm(1,1) != 17 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 18  0  0 )\n( 15 17 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != 12 || lower(1,1) != 18 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 15 || lower(2,1) != 17 || lower(2,2) != 11 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( 12 18  0  0 )\n"
                                     "( 15 17 11  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 11 12  0 )
   // (  7  0  3  0 )      (  7 13 14  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major dense matrix assignment test 3";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse matrix assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 18  0  0 )
   // (  7  0  3  0 )      (  7 14 11  0 )
   // ( -2  0  1  5 )      ( -2 15 19  5 )
   {
      test_ = "Row-major sparse matrix assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 6UL );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;
      mat.insert( 0UL, 1UL, 0 );

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) != 14 || sm(2,1) != 11 ||
          sm(3,0) != 15 || sm(3,1) != 19 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n( 14 11 )\n( 15 19 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 18 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 14 || lower(2,2) != 11 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 15 || lower(3,2) != 19 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 18  0  0 )\n"
                                     "(  7 14 11  0 )\n"
                                     "( -2 15 19  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( 12 18  0  0 )
   // (  7  0  3  0 )      ( 15 17 11  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major sparse matrix assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 6UL );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;
      mat.insert( 0UL, 3UL, 0 );

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 15 || sm(1,1) != 17 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 18  0  0 )\n( 15 17 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != 12 || lower(1,1) != 18 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 15 || lower(2,1) != 17 || lower(2,2) != 11 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( 12 18  0  0 )\n"
                                     "( 15 17 11  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 11 12  0 )
   // (  7  0  3  0 )      (  7 13 14  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major sparse matrix assignment test 3";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 18  0  0 )
   // (  7  0  3  0 )      (  7 14 11  0 )
   // ( -2  0  1  5 )      ( -2 15 19  5 )
   {
      test_ = "Column-major sparse matrix assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 6UL );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;
      mat.insert( 0UL, 1UL, 0 );

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) != 14 || sm(2,1) != 11 ||
          sm(3,0) != 15 || sm(3,1) != 19 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n( 14 11 )\n( 15 19 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 18 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 14 || lower(2,2) != 11 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 15 || lower(3,2) != 19 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 18  0  0 )\n"
                                     "(  7 14 11  0 )\n"
                                     "( -2 15 19  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( 12 18  0  0 )
   // (  7  0  3  0 )      ( 15 17 11  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major sparse matrix assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 6UL );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;
      mat.insert( 0UL, 3UL, 0 );

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 15 || sm(1,1) != 17 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 18  0  0 )\n( 15 17 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != 12 || lower(1,1) != 18 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 15 || lower(2,1) != 17 || lower(2,2) != 11 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( 12 18  0  0 )\n"
                                     "( 15 17 11  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 11 12  0 )
   // (  7  0  3  0 )      (  7 13 14  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major sparse matrix assignment test 3";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the addition assignment to a submatrix of a LowerMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment to a submatrix of a LowerMatrix. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename LT >  // Type of the lower matrix
void SubmatrixTest::testAddAssign()
{
   //=====================================================================================
   // Dense matrix addition assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 20  0  0 )
   // (  7  0  3  0 )      (  7 14 14  0 )
   // ( -2  0  1  5 )      ( -2 15 20  5 )
   {
      test_ = "Row-major dense matrix addition assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 20 || sm(1,1) !=  0 ||
          sm(2,0) != 14 || sm(2,1) != 14 ||
          sm(3,0) != 15 || sm(3,1) != 20 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 20  0 )\n( 14 14 )\n( 15 20 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 20 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 14 || lower(2,2) != 14 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 15 || lower(3,2) != 20 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 20  0  0 )\n"
                                     "(  7 14 14  0 )\n"
                                     "( -2 15 20  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  (  8 20  0  0 )
   // (  7  0  3  0 )      ( 22 17 14  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major dense matrix addition assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) !=  8 || sm(0,1) != 20 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 22 || sm(1,1) != 17 || sm(1,2) != 14 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  8 20  0  0 )\n( 22 17 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) !=  8 || lower(1,1) != 20 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 22 || lower(2,1) != 17 || lower(2,2) != 14 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "(  8 20  0  0 )\n"
                                     "( 22 17 14  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 13 12  0 )
   // (  7  0  3  0 )      (  7 13 17  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major dense matrix addition assignment test 3";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 20  0  0 )
   // (  7  0  3  0 )      (  7 14 14  0 )
   // ( -2  0  1  5 )      ( -2 15 20  5 )
   {
      test_ = "Column-major dense matrix addition assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 20 || sm(1,1) !=  0 ||
          sm(2,0) != 14 || sm(2,1) != 14 ||
          sm(3,0) != 15 || sm(3,1) != 20 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 20  0 )\n( 14 14 )\n( 15 20 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 20 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 14 || lower(2,2) != 14 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 15 || lower(3,2) != 20 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 20  0  0 )\n"
                                     "(  7 14 14  0 )\n"
                                     "( -2 15 20  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  (  8 20  0  0 )
   // (  7  0  3  0 )      ( 22 17 14  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major dense matrix addition assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) !=  8 || sm(0,1) != 20 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 22 || sm(1,1) != 17 || sm(1,2) != 14 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 18  0  0 )\n( 15 17 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) !=  8 || lower(1,1) != 20 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 22 || lower(2,1) != 17 || lower(2,2) != 14 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "(  8 20  0  0 )\n"
                                     "( 22 17 14  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 13 12  0 )
   // (  7  0  3  0 )      (  7 13 17  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major dense matrix addition assignment test 3";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 20  0  0 )
   // (  7  0  3  0 )      (  7 14 14  0 )
   // ( -2  0  1  5 )      ( -2 15 20  5 )
   {
      test_ = "Row-major sparse matrix addition assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 5UL );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 20 || sm(1,1) !=  0 ||
          sm(2,0) != 14 || sm(2,1) != 14 ||
          sm(3,0) != 15 || sm(3,1) != 20 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 20  0 )\n( 14 14 )\n( 15 20 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 20 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 14 || lower(2,2) != 14 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 15 || lower(3,2) != 20 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 20  0  0 )\n"
                                     "(  7 14 14  0 )\n"
                                     "( -2 15 20  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  (  8 20  0  0 )
   // (  7  0  3  0 )      ( 22 17 14  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major sparse matrix addition assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 5UL );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) !=  8 || sm(0,1) != 20 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 22 || sm(1,1) != 17 || sm(1,2) != 14 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  8 20  0  0 )\n( 22 17 14  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) !=  8 || lower(1,1) != 20 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 22 || lower(2,1) != 17 || lower(2,2) != 14 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "(  8 20  0  0 )\n"
                                     "( 22 17 14  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 13 12  0 )
   // (  7  0  3  0 )      (  7 13 17  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major sparse matrix addition assignment test 3";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 20  0  0 )
   // (  7  0  3  0 )      (  7 14 14  0 )
   // ( -2  0  1  5 )      ( -2 15 20  5 )
   {
      test_ = "Column-major sparse matrix addition assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 5UL );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 20 || sm(1,1) !=  0 ||
          sm(2,0) != 14 || sm(2,1) != 14 ||
          sm(3,0) != 15 || sm(3,1) != 20 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 20  0 )\n( 14 14 )\n( 15 20 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 20 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != 14 || lower(2,2) != 14 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != 15 || lower(3,2) != 20 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 20  0  0 )\n"
                                     "(  7 14 14  0 )\n"
                                     "( -2 15 20  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  (  8 20  0  0 )
   // (  7  0  3  0 )      ( 22 17 14  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major sparse matrix addition assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 5UL );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) !=  8 || sm(0,1) != 20 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 22 || sm(1,1) != 17 || sm(1,2) != 14 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 18  0  0 )\n( 15 17 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) !=  8 || lower(1,1) != 20 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 22 || lower(2,1) != 17 || lower(2,2) != 14 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "(  8 20  0  0 )\n"
                                     "( 22 17 14  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 13 12  0 )
   // (  7  0  3  0 )      (  7 13 17  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major sparse matrix addition assignment test 3";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the subtraction assignment to a submatrix of a LowerMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment to a submatrix of a LowerMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename LT >  // Type of the lower matrix
void SubmatrixTest::testSubAssign()
{
   //=====================================================================================
   // Dense matrix subtraction assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1   0   0   0 )
   // ( -4  2  0  0 )  =>  ( -4 -16   0   0 )
   // (  7  0  3  0 )      (  7 -14  -8   0 )
   // ( -2  0  1  5 )      ( -2 -15 -18   5 )
   {
      test_ = "Row-major dense matrix subtraction assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=   0 || sm(0,1) !=   0 ||
          sm(1,0) != -16 || sm(1,1) !=   0 ||
          sm(2,0) != -14 || sm(2,1) !=  -8 ||
          sm(3,0) != -15 || sm(3,1) != -18 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(   0   0 )\n( -16   0 )\n( -14  -8 )\n( -15 -18 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=   0 || lower(0,2) !=   0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -16 || lower(1,2) !=   0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -14 || lower(2,2) !=  -8 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != -15 || lower(3,2) != -18 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1   0   0   0 )\n"
                                     "( -4 -16   0   0 )\n"
                                     "(  7 -14  -8   0 )\n"
                                     "( -2 -15 -18   5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (   1   0   0   0 )
   // ( -4  2  0  0 )  =>  ( -16 -16   0   0 )
   // (  7  0  3  0 )      (  -8 -17  -8   0 )
   // ( -2  0  1  5 )      (  -2   0   1   5 )
   {
      test_ = "Row-major dense matrix subtraction assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) != -16 || sm(0,1) != -16 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) !=  -8 || sm(1,1) != -17 || sm(1,2) != -8 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -16 -16  0  0 )\n(  -8 -17 -8  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=   1 || lower(0,1) !=   0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -16 || lower(1,1) != -16 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  -8 || lower(2,1) != -17 || lower(2,2) != -8 || lower(2,3) != 0 ||
          lower(3,0) !=  -2 || lower(3,1) !=   0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(   1   0  0  0 )\n"
                                     "( -16 -16  0  0 )\n"
                                     "(  -8 -17 -8  0 )\n"
                                     "(  -2  0   1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1   0   0  0 )
   // ( -4  2  0  0 )  =>  ( -4  -9 -12  0 )
   // (  7  0  3  0 )      (  7 -13 -11  0 )
   // ( -2  0  1  5 )      ( -2   0   1  5 )
   {
      test_ = "Row-major dense matrix subtraction assignment test 3";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0  0 )      (  1   0   0   0 )
   // ( -4  2  0  0 )  =>  ( -4 -16   0   0 )
   // (  7  0  3  0 )      (  7 -14  -8   0 )
   // ( -2  0  1  5 )      ( -2 -15 -18   5 )
   {
      test_ = "Column-major dense matrix subtraction assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=   0 || sm(0,1) !=   0 ||
          sm(1,0) != -16 || sm(1,1) !=   0 ||
          sm(2,0) != -14 || sm(2,1) !=  -8 ||
          sm(3,0) != -15 || sm(3,1) != -18 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(   0   0 )\n( -16   0 )\n( -14  -8 )\n( -15 -18 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=   0 || lower(0,2) !=   0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -16 || lower(1,2) !=   0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -14 || lower(2,2) !=  -8 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != -15 || lower(3,2) != -18 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1   0   0   0 )\n"
                                     "( -4 -16   0   0 )\n"
                                     "(  7 -14  -8   0 )\n"
                                     "( -2 -15 -18   5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (   1   0   0   0 )
   // ( -4  2  0  0 )  =>  ( -16 -16   0   0 )
   // (  7  0  3  0 )      (  -8 -17  -8   0 )
   // ( -2  0  1  5 )      (  -2   0   1   5 )
   {
      test_ = "Column-major dense matrix subtraction assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) != -16 || sm(0,1) != -16 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) !=  -8 || sm(1,1) != -17 || sm(1,2) != -8 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -16 -16  0  0 )\n(  -8 -17 -8  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=   1 || lower(0,1) !=   0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -16 || lower(1,1) != -16 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  -8 || lower(2,1) != -17 || lower(2,2) != -8 || lower(2,3) != 0 ||
          lower(3,0) !=  -2 || lower(3,1) !=   0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(   1   0  0  0 )\n"
                                     "( -16 -16  0  0 )\n"
                                     "(  -8 -17 -8  0 )\n"
                                     "(  -2  0   1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1   0   0  0 )
   // ( -4  2  0  0 )  =>  ( -4  -9 -12  0 )
   // (  7  0  3  0 )      (  7 -13 -11  0 )
   // ( -2  0  1  5 )      ( -2   0   1  5 )
   {
      test_ = "Column-major dense matrix subtraction assignment test 3";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse matrix subtraction assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1   0   0   0 )
   // ( -4  2  0  0 )  =>  ( -4 -16   0   0 )
   // (  7  0  3  0 )      (  7 -14  -8   0 )
   // ( -2  0  1  5 )      ( -2 -15 -18   5 )
   {
      test_ = "Row-major sparse matrix subtraction assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 5UL );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=   0 || sm(0,1) !=   0 ||
          sm(1,0) != -16 || sm(1,1) !=   0 ||
          sm(2,0) != -14 || sm(2,1) !=  -8 ||
          sm(3,0) != -15 || sm(3,1) != -18 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(   0   0 )\n( -16   0 )\n( -14  -8 )\n( -15 -18 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=   0 || lower(0,2) !=   0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -16 || lower(1,2) !=   0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -14 || lower(2,2) !=  -8 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != -15 || lower(3,2) != -18 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1   0   0   0 )\n"
                                     "( -4 -16   0   0 )\n"
                                     "(  7 -14  -8   0 )\n"
                                     "( -2 -15 -18   5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (   1   0   0   0 )
   // ( -4  2  0  0 )  =>  ( -16 -16   0   0 )
   // (  7  0  3  0 )      (  -8 -17  -8   0 )
   // ( -2  0  1  5 )      (  -2   0   1   5 )
   {
      test_ = "Row-major sparse matrix subtraction assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 5UL );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) != -16 || sm(0,1) != -16 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) !=  -8 || sm(1,1) != -17 || sm(1,2) != -8 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -16 -16  0  0 )\n(  -8 -17 -8  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=   1 || lower(0,1) !=   0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -16 || lower(1,1) != -16 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  -8 || lower(2,1) != -17 || lower(2,2) != -8 || lower(2,3) != 0 ||
          lower(3,0) !=  -2 || lower(3,1) !=   0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(   1   0  0  0 )\n"
                                     "( -16 -16  0  0 )\n"
                                     "(  -8 -17 -8  0 )\n"
                                     "(  -2  0   1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1   0   0  0 )
   // ( -4  2  0  0 )  =>  ( -4  -9 -12  0 )
   // (  7  0  3  0 )      (  7 -13 -11  0 )
   // ( -2  0  1  5 )      ( -2   0   1  5 )
   {
      test_ = "Row-major sparse matrix subtraction assignment test 3";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // (  1  0  0  0 )      (  1   0   0   0 )
   // ( -4  2  0  0 )  =>  ( -4 -16   0   0 )
   // (  7  0  3  0 )      (  7 -14  -8   0 )
   // ( -2  0  1  5 )      ( -2 -15 -18   5 )
   {
      test_ = "Column-major sparse matrix subtraction assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 5UL );
      mat(1,0) = 18;
      mat(2,0) = 14;
      mat(2,1) = 11;
      mat(3,0) = 15;
      mat(3,1) = 19;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( lower,  4UL );
      checkColumns ( lower,  4UL );
      checkNonZeros( lower, 10UL );

      if( sm(0,0) !=   0 || sm(0,1) !=   0 ||
          sm(1,0) != -16 || sm(1,1) !=   0 ||
          sm(2,0) != -14 || sm(2,1) !=  -8 ||
          sm(3,0) != -15 || sm(3,1) != -18 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(   0   0 )\n( -16   0 )\n( -14  -8 )\n( -15 -18 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=   0 || lower(0,2) !=   0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != -16 || lower(1,2) !=   0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) != -14 || lower(2,2) !=  -8 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) != -15 || lower(3,2) != -18 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1   0   0   0 )\n"
                                     "( -4 -16   0   0 )\n"
                                     "(  7 -14  -8   0 )\n"
                                     "( -2 -15 -18   5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (   1   0   0   0 )
   // ( -4  2  0  0 )  =>  ( -16 -16   0   0 )
   // (  7  0  3  0 )      (  -8 -17  -8   0 )
   // ( -2  0  1  5 )      (  -2   0   1   5 )
   {
      test_ = "Column-major sparse matrix subtraction assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 5UL );
      mat(0,0) = 12;
      mat(0,1) = 18;
      mat(1,0) = 15;
      mat(1,1) = 17;
      mat(1,2) = 11;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 9UL );

      if( sm(0,0) != -16 || sm(0,1) != -16 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) !=  -8 || sm(1,1) != -17 || sm(1,2) != -8 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -16 -16  0  0 )\n(  -8 -17 -8  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=   1 || lower(0,1) !=   0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -16 || lower(1,1) != -16 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  -8 || lower(2,1) != -17 || lower(2,2) != -8 || lower(2,3) != 0 ||
          lower(3,0) !=  -2 || lower(3,1) !=   0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(   1   0  0  0 )\n"
                                     "( -16 -16  0  0 )\n"
                                     "(  -8 -17 -8  0 )\n"
                                     "(  -2  0   1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1   0   0  0 )
   // ( -4  2  0  0 )  =>  ( -4  -9 -12  0 )
   // (  7  0  3  0 )      (  7 -13 -11  0 )
   // ( -2  0  1  5 )      ( -2   0   1  5 )
   {
      test_ = "Column-major sparse matrix subtraction assignment test 3";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) = 14;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Schur product assignment to a submatrix of a LowerMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment to a submatrix of a LowerMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename LT >  // Type of the lower matrix
void SubmatrixTest::testSchurAssign()
{
   //=====================================================================================
   // Dense matrix Schur product assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 12  0  0 )
   // (  7  0  3  0 )      (  7  0 12  0 )
   // ( -2  0  1  5 )      ( -2  0  0  5 )
   {
      test_ = "Row-major dense matrix Schur product assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(0,1) = 9;
      mat(1,0) = 6;
      mat(2,1) = 4;
      mat(3,0) = 9;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 7UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 12 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 12 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 12  0 )\n(  0 12 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  0 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 12  0  0 )\n"
                                     "(  7  0 12  0 )\n"
                                     "( -2  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  (  0 12  0  0 )
   // (  7  0  3  0 )      ( 21  0 12  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major dense matrix Schur product assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 6;
      mat(0,3) = 9;
      mat(1,0) = 3;
      mat(1,2) = 4;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 7UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 21 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0  0 )\n( 21  0 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) !=  0 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 21 || lower(2,1) !=  0 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 21  0 12  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 12  0  0 )
   // (  7  0  3  0 )      (  7  0 12  0 )
   // ( -2  0  1  5 )      ( -2  0  0  5 )
   {
      test_ = "Column-major dense matrix Schur product assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(0,1) = 9;
      mat(1,0) = 6;
      mat(2,1) = 4;
      mat(3,0) = 9;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 7UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 12 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 12 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 12  0 )\n(  0 12 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  0 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 12  0  0 )\n"
                                     "(  7  0 12  0 )\n"
                                     "( -2  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  (  0 12  0  0 )
   // (  7  0  3  0 )      ( 21  0 12  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major dense matrix Schur product assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 6;
      mat(0,3) = 9;
      mat(1,0) = 3;
      mat(1,2) = 4;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 7UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 21 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0  0 )\n( 21  0 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) !=  0 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 21 || lower(2,1) !=  0 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 21  0 12  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix Schur product assignment
   //=====================================================================================

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 12  0  0 )
   // (  7  0  3  0 )      (  7  0 12  0 )
   // ( -2  0  1  5 )      ( -2  0  0  5 )
   {
      test_ = "Row-major sparse matrix Schur product assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 4UL );
      mat(0,1) = 9;
      mat(1,0) = 6;
      mat(2,1) = 4;
      mat(3,0) = 9;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 7UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 12 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 12 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 12  0 )\n(  0 12 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  0 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 12  0  0 )\n"
                                     "(  7  0 12  0 )\n"
                                     "( -2  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  (  0 12  0  0 )
   // (  7  0  3  0 )      ( 21  0 12  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Row-major sparse matrix Schur product assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 6;
      mat(0,3) = 9;
      mat(1,0) = 3;
      mat(1,2) = 4;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 7UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 21 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0  0 )\n( 21  0 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) !=  0 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 21 || lower(2,1) !=  0 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 21  0 12  0 )\n"
                                     "( -2  0  1  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  ( -4 12  0  0 )
   // (  7  0  3  0 )      (  7  0 12  0 )
   // ( -2  0  1  5 )      ( -2  0  0  5 )
   {
      test_ = "Column-major sparse matrix Schur product assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 4UL );
      mat(0,1) = 9;
      mat(1,0) = 6;
      mat(2,1) = 4;
      mat(3,0) = 9;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 7UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 12 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 12 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 12  0 )\n(  0 12 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) != -4 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) !=  7 || lower(2,1) !=  0 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  0 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "( -4 12  0  0 )\n"
                                     "(  7  0 12  0 )\n"
                                     "( -2  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // (  1  0  0  0 )      (  1  0  0  0 )
   // ( -4  2  0  0 )  =>  (  0 12  0  0 )
   // (  7  0  3  0 )      ( 21  0 12  0 )
   // ( -2  0  1  5 )      ( -2  0  1  5 )
   {
      test_ = "Column-major sparse matrix Schur product assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 6;
      mat(0,3) = 9;
      mat(1,0) = 3;
      mat(1,2) = 4;

      LT lower;
      init( lower );

      auto sm = submatrix( lower, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( lower, 4UL );
      checkColumns ( lower, 4UL );
      checkNonZeros( lower, 7UL );

      if( sm(0,0) !=  0 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 21 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0 12  0  0 )\n( 21  0 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( lower(0,0) !=  1 || lower(0,1) !=  0 || lower(0,2) !=  0 || lower(0,3) != 0 ||
          lower(1,0) !=  0 || lower(1,1) != 12 || lower(1,2) !=  0 || lower(1,3) != 0 ||
          lower(2,0) != 21 || lower(2,1) !=  0 || lower(2,2) != 12 || lower(2,3) != 0 ||
          lower(3,0) != -2 || lower(3,1) !=  0 || lower(3,2) !=  1 || lower(3,3) != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << lower << "\n"
             << "   Expected result:\n(  1  0  0  0 )\n"
                                     "(  0 12  0  0 )\n"
                                     "( 21  0 12  0 )\n"
                                     "( -2  0  1  5 )\n";
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
void SubmatrixTest::checkRows( const Type& matrix, size_t expectedRows ) const
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
void SubmatrixTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
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
void SubmatrixTest::checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const
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
/*!\brief Initializing the given lower triangular matrix.
//
// \return void
//
// This function is called before each test case to initialize the given lower triangular matrix.
*/
template< typename LT >
void SubmatrixTest::init( LT& lower )
{
   lower.resize( 4UL );
   lower(0,0) =  1;
   lower(1,0) = -4;
   lower(1,1) =  2;
   lower(2,0) =  7;
   lower(2,1) =  0;
   lower(2,2) =  3;
   lower(3,0) = -2;
   lower(3,1) =  0;
   lower(3,2) =  1;
   lower(3,3) =  5;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the assignment to a submatrix of a LowerMatrix.
//
// \return void
*/
void runTest()
{
   SubmatrixTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the LowerMatrix submatrix test.
*/
#define RUN_LOWERMATRIX_SUBMATRIX_TEST \
   blazetest::mathtest::adaptors::lowermatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace lowermatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
