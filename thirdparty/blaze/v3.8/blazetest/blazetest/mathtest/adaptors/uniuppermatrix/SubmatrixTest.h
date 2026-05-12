//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/uniuppermatrix/SubmatrixTest.h
//  \brief Header file for the UniUpperMatrix submatrix test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_UNIUPPERMATRIX_SUBMATRIXTEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_UNIUPPERMATRIX_SUBMATRIXTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Submatrix.h>
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
/*!\brief Auxiliary class for assignment tests to a submatrix of a UniUpperMatrix.
//
// This class performs assignment tests to a submatrix of a UniUpperMatrix. It performs a series
// of both compile time as well as runtime tests.
*/
class SubmatrixTest
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
   template< typename UT > void testAssignment ();
   template< typename UT > void testAddAssign  ();
   template< typename UT > void testSubAssign  ();
   template< typename UT > void testSchurAssign();

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
/*!\brief Test of the assignment to a submatrix of a UniUpperMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to a submatrix of a UniUpperMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename UT >  // Type of the uniupper matrix
void SubmatrixTest::testAssignment()
{
   //=====================================================================================
   // Dense matrix assignment
   //=====================================================================================

   // ( 1 -4  7 -2 )      ( 1 12 15 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 17  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major dense matrix assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,0) =  1;
      mat(1,1) = 17;
      mat(2,1) =  1;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != 12 || sm(0,1) != 15 ||
          sm(1,0) !=  1 || sm(1,1) != 17 ||
          sm(2,0) !=  0 || sm(2,1) !=  1 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 15 )\n(  1 17 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 12 || upper(0,2) != 15 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 12 15 -2 )\n"
                                     "( 0  1 17  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 14 10 )
   // ( 0  0  1  1 )      ( 0  0  1 16 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major dense matrix assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) =  1;
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,2) =  1;
      mat(1,3) = 16;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 14 || sm(0,3) != 10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  1 || sm(1,3) != 16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 14 10 )\n( 0  0  1 16 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 14 || upper(1,3) != 10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) != 16 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1 14 10 )\n"
                                     "( 0  0  1 16 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0  0  0  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Row-major dense matrix assignment test 3";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;
      mat(1,0) =  0;
      mat(1,1) =  0;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0 13  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Row-major dense matrix assignment test 4";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) =  1;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 12 15 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 17  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major dense matrix assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,0) =  1;
      mat(1,1) = 17;
      mat(2,1) =  1;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != 12 || sm(0,1) != 15 ||
          sm(1,0) !=  1 || sm(1,1) != 17 ||
          sm(2,0) !=  0 || sm(2,1) !=  1 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 15 )\n(  1 17 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 12 || upper(0,2) != 15 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 12 15 -2 )\n"
                                     "( 0  1 17  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 14 10 )
   // ( 0  0  1  1 )      ( 0  0  1 16 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major dense matrix assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) =  1;
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,2) =  1;
      mat(1,3) = 16;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 14 || sm(0,3) != 10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  1 || sm(1,3) != 16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 14 10 )\n( 0  0  1 16 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 14 || upper(1,3) != 10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) != 16 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1 14 10 )\n"
                                     "( 0  0  1 16 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0  0  0  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Column-major dense matrix assignment test 3";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;
      mat(1,0) =  0;
      mat(1,1) =  0;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0 13  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Column-major dense matrix assignment test 4";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) =  1;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse matrix assignment
   //=====================================================================================

   // ( 1 -4  7 -2 )      ( 1 12 15 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 17  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major sparse matrix assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 6UL );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,0) =  1;
      mat(1,1) = 17;
      mat(2,1) =  1;
      mat.insert( 3UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != 12 || sm(0,1) != 15 ||
          sm(1,0) !=  1 || sm(1,1) != 17 ||
          sm(2,0) !=  0 || sm(2,1) !=  1 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 15 )\n(  1 17 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 12 || upper(0,2) != 15 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 12 15 -2 )\n"
                                     "( 0  1 17  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 14 10 )
   // ( 0  0  1  1 )      ( 0  0  1 16 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major sparse matrix assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 6UL );
      mat(0,1) =  1;
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,2) =  1;
      mat(1,3) = 16;
      mat.insert( 1UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 14 || sm(0,3) != 10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  1 || sm(1,3) != 16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 14 10 )\n( 0  0  1 16 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 14 || upper(1,3) != 10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) != 16 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1 14 10 )\n"
                                     "( 0  0  1 16 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0  0  0  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Row-major sparse matrix assignment test 3";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0 13  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Row-major sparse matrix assignment test 4";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  1;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) =  1;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 12 15 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 17  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major sparse matrix assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 6UL );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,0) =  1;
      mat(1,1) = 17;
      mat(2,1) =  1;
      mat.insert( 3UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != 12 || sm(0,1) != 15 ||
          sm(1,0) !=  1 || sm(1,1) != 17 ||
          sm(2,0) !=  0 || sm(2,1) !=  1 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 12 15 )\n(  1 17 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 12 || upper(0,2) != 15 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 12 15 -2 )\n"
                                     "( 0  1 17  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 14 10 )
   // ( 0  0  1  1 )      ( 0  0  1 16 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major sparse matrix assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 6UL );
      mat(0,1) =  1;
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,2) =  1;
      mat(1,3) = 16;
      mat.insert( 1UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 14 || sm(0,3) != 10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  1 || sm(1,3) != 16 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 14 10 )\n( 0  0  1 16 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 14 || upper(1,3) != 10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) != 16 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1 14 10 )\n"
                                     "( 0  0  1 16 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0  0  0  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Column-major sparse matrix assignment test 3";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0 13  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Column-major sparse matrix assignment test 4";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  1;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) =  1;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the addition assignment to a submatrix of a UniUpperMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment to a submatrix of a UniUpperMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename UT >  // Type of the uniupper matrix
void SubmatrixTest::testAddAssign()
{
   //=====================================================================================
   // Dense matrix addition assignment
   //=====================================================================================

   // ( 1 -4  7 -2 )      ( 1  8 22 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 17  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major dense matrix addition assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,1) = 17;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != 8 || sm(0,1) != 22 ||
          sm(1,0) != 1 || sm(1,1) != 17 ||
          sm(2,0) != 0 || sm(2,1) !=  1 ||
          sm(3,0) != 0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  8 22 )\n(  1 17 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 8 || upper(0,2) != 22 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1  8 22 -2 )\n"
                                     "( 0  1 17  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 14 10 )
   // ( 0  0  1  1 )      ( 0  0  1 17 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major dense matrix addition assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,3) = 16;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 14 || sm(0,3) != 10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  1 || sm(1,3) != 17 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 14 10 )\n( 0  0  1 17 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 14 || upper(1,3) != 10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) != 17 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1 14 10 )\n"
                                     "( 0  0  1 17 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  2 12  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Row-major dense matrix addition assignment test 3";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;
      mat(1,0) =  0;
      mat(1,1) =  0;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0 13  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Row-major dense matrix addition assignment test 4";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) =  0;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) =  0;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1  8 22 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 17  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major dense matrix addition assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,1) = 17;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != 8 || sm(0,1) != 22 ||
          sm(1,0) != 1 || sm(1,1) != 17 ||
          sm(2,0) != 0 || sm(2,1) !=  1 ||
          sm(3,0) != 0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  8 22 )\n(  1 17 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 8 || upper(0,2) != 22 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1  8 22 -2 )\n"
                                     "( 0  1 17  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 14 10 )
   // ( 0  0  1  1 )      ( 0  0  1 17 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major dense matrix addition assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,3) = 16;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 14 || sm(0,3) != 10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  1 || sm(1,3) != 17 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 14 10 )\n( 0  0  1 17 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 14 || upper(1,3) != 10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) != 17 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1 14 10 )\n"
                                     "( 0  0  1 17 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  2 12  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Column-major dense matrix addition assignment test 3";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;
      mat(1,0) =  0;
      mat(1,1) =  0;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0 13  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Column-major dense matrix addition assignment test 4";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) =  0;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) =  0;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse matrix addition assignment
   //=====================================================================================

   // ( 1 -4  7 -2 )      ( 1  8 22 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 17  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major sparse matrix addition assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 4UL );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,1) = 17;
      mat.insert( 3UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != 8 || sm(0,1) != 22 ||
          sm(1,0) != 1 || sm(1,1) != 17 ||
          sm(2,0) != 0 || sm(2,1) !=  1 ||
          sm(3,0) != 0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  8 22 )\n(  1 17 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 8 || upper(0,2) != 22 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1  8 22 -2 )\n"
                                     "( 0  1 17  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 14 10 )
   // ( 0  0  1  1 )      ( 0  0  1 17 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major sparse matrix addition assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 4UL );
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,3) = 16;
      mat.insert( 1UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 14 || sm(0,3) != 10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  1 || sm(1,3) != 17 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 14 10 )\n( 0  0  1 17 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 14 || upper(1,3) != 10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) != 17 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1 14 10 )\n"
                                     "( 0  0  1 17 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  2 12  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Row-major sparse matrix addition assignment test 3";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0 13  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Row-major sparse matrix addition assignment test 4";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 2UL );
      mat(0,1) = 12;
      mat(1,0) = 13;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1  8 22 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 17  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major sparse matrix addition assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 4UL );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,1) = 17;
      mat.insert( 3UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != 8 || sm(0,1) != 22 ||
          sm(1,0) != 1 || sm(1,1) != 17 ||
          sm(2,0) != 0 || sm(2,1) !=  1 ||
          sm(3,0) != 0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  8 22 )\n(  1 17 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 8 || upper(0,2) != 22 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) != 17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1  8 22 -2 )\n"
                                     "( 0  1 17  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 14 10 )
   // ( 0  0  1  1 )      ( 0  0  1 17 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major sparse matrix addition assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 4UL );
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,3) = 16;
      mat.insert( 1UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 14 || sm(0,3) != 10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=  1 || sm(1,3) != 17 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 14 10 )\n( 0  0  1 17 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=  7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 14 || upper(1,3) != 10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=  1 || upper(2,3) != 17 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1 14 10 )\n"
                                     "( 0  0  1 17 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  2 12  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Column-major sparse matrix addition assignment test 3";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 12  0 )
   // ( 0  0  1  1 )      ( 0 13  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Column-major sparse matrix addition assignment test 4";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 2UL );
      mat(0,1) = 12;
      mat(1,0) = 13;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the subtraction assignment to a submatrix of a UniUpperMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment to a submatrix of a UniUpperMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename UT >  // Type of the uniupper matrix
void SubmatrixTest::testSubAssign()
{
   //=====================================================================================
   // Dense matrix subtraction assignment
   //=====================================================================================

   // ( 1 -4  7 -2 )      ( 1 -16  -8 -2 )
   // ( 0  1  0  0 )  =>  ( 0   1 -17  0 )
   // ( 0  0  1  1 )      ( 0   0   1  1 )
   // ( 0  0  0  1 )      ( 0   0   0  1 )
   {
      test_ = "Row-major dense matrix subtraction assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,1) = 17;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != -16 || sm(0,1) !=  -8 ||
          sm(1,0) !=   1 || sm(1,1) != -17 ||
          sm(2,0) !=   0 || sm(2,1) !=   1 ||
          sm(3,0) !=   0 || sm(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -16  -8 )\n(   1 -17 )\n(   0   1 )\n(   0   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -16 || upper(0,2) !=  -8 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=   1 || upper(1,2) != -17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=   0 || upper(2,2) !=   1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=   0 || upper(3,2) !=   0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -16  -8 -2 )\n"
                                     "( 0   1 -17  0 )\n"
                                     "( 0   0   1  1 )\n"
                                     "( 0   0   0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4   7  -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 -14 -10 )
   // ( 0  0  1  1 )      ( 0  0   1 -15 )
   // ( 0  0  0  1 )      ( 0  0   0   1 )
   {
      test_ = "Row-major dense matrix subtraction assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,3) = 16;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != -14 || sm(0,3) != -10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=   1 || sm(1,3) != -15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 -14 -10 )\n( 0  0   1 -15 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=   7 || upper(0,3) !=  -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -14 || upper(1,3) != -10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=   1 || upper(2,3) != -15 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=   0 || upper(3,3) !=   1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4   7  -2 )\n"
                                     "( 0  1 -14 -10 )\n"
                                     "( 0  0   1 -15 )\n"
                                     "( 0  0   0   1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4   7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  0 -12  0 )
   // ( 0  0  1  1 )      ( 0  0   1  1 )
   // ( 0  0  0  1 )      ( 0  0   0  5 )
   {
      test_ = "Row-major dense matrix subtraction assignment test 3";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;
      mat(1,0) =  0;
      mat(1,1) =  0;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1  -4   7 -2 )
   // ( 0  1  0  0 )  =>  ( 0   1 -12  0 )
   // ( 0  0  1  1 )      ( 0 -13   1  1 )
   // ( 0  0  0  1 )      ( 0   0   0  5 )
   {
      test_ = "Row-major dense matrix subtraction assignment test 4";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) =  0;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) =  0;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 -16  -8 -2 )
   // ( 0  1  0  0 )  =>  ( 0   1 -17  0 )
   // ( 0  0  1  1 )      ( 0   0   1  1 )
   // ( 0  0  0  1 )      ( 0   0   0  1 )
   {
      test_ = "Column-major dense matrix subtraction assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,1) = 17;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != -16 || sm(0,1) !=  -8 ||
          sm(1,0) !=   1 || sm(1,1) != -17 ||
          sm(2,0) !=   0 || sm(2,1) !=   1 ||
          sm(3,0) !=   0 || sm(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -16  -8 )\n(   1 -17 )\n(   0   1 )\n(   0   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -16 || upper(0,2) !=  -8 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=   1 || upper(1,2) != -17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=   0 || upper(2,2) !=   1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=   0 || upper(3,2) !=   0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -16  -8 -2 )\n"
                                     "( 0   1 -17  0 )\n"
                                     "( 0   0   1  1 )\n"
                                     "( 0   0   0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4   7  -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 -14 -10 )
   // ( 0  0  1  1 )      ( 0  0   1 -15 )
   // ( 0  0  0  1 )      ( 0  0   0   1 )
   {
      test_ = "Column-major dense matrix subtraction assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,3) = 16;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != -14 || sm(0,3) != -10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=   1 || sm(1,3) != -15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 -14 -10 )\n( 0  0   1 -15 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=   7 || upper(0,3) !=  -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -14 || upper(1,3) != -10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=   1 || upper(2,3) != -15 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=   0 || upper(3,3) !=   1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4   7  -2 )\n"
                                     "( 0  1 -14 -10 )\n"
                                     "( 0  0   1 -15 )\n"
                                     "( 0  0   0   1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4   7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  0 -12  0 )
   // ( 0  0  1  1 )      ( 0  0   1  1 )
   // ( 0  0  0  1 )      ( 0  0   0  5 )
   {
      test_ = "Column-major dense matrix subtraction assignment test 3";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;
      mat(1,0) =  0;
      mat(1,1) =  0;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1  -4   7 -2 )
   // ( 0  1  0  0 )  =>  ( 0   1 -12  0 )
   // ( 0  0  1  1 )      ( 0 -13   1  1 )
   // ( 0  0  0  1 )      ( 0   0   0  5 )
   {
      test_ = "Column-major dense matrix subtraction assignment test 4";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) =  0;
      mat(0,1) = 12;
      mat(1,0) = 13;
      mat(1,1) =  0;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse matrix subtraction assignment
   //=====================================================================================

   // ( 1 -4  7 -2 )      ( 1 -16  -8 -2 )
   // ( 0  1  0  0 )  =>  ( 0   1 -17  0 )
   // ( 0  0  1  1 )      ( 0   0   1  1 )
   // ( 0  0  0  1 )      ( 0   0   0  1 )
   {
      test_ = "Row-major sparse matrix subtraction assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 4UL );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,1) = 17;
      mat.insert( 3UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != -16 || sm(0,1) !=  -8 ||
          sm(1,0) !=   1 || sm(1,1) != -17 ||
          sm(2,0) !=   0 || sm(2,1) !=   1 ||
          sm(3,0) !=   0 || sm(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -16  -8 )\n(   1 -17 )\n(   0   1 )\n(   0   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -16 || upper(0,2) !=  -8 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=   1 || upper(1,2) != -17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=   0 || upper(2,2) !=   1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=   0 || upper(3,2) !=   0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -16  -8 -2 )\n"
                                     "( 0   1 -17  0 )\n"
                                     "( 0   0   1  1 )\n"
                                     "( 0   0   0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4   7  -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 -14 -10 )
   // ( 0  0  1  1 )      ( 0  0   1 -15 )
   // ( 0  0  0  1 )      ( 0  0   0   1 )
   {
      test_ = "Row-major sparse matrix subtraction assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 4UL );
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,3) = 16;
      mat.insert( 1UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != -14 || sm(0,3) != -10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=   1 || sm(1,3) != -15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 -14 -10 )\n( 0  0   1 -15 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=   7 || upper(0,3) !=  -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -14 || upper(1,3) != -10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=   1 || upper(2,3) != -15 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=   0 || upper(3,3) !=   1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4   7  -2 )\n"
                                     "( 0  1 -14 -10 )\n"
                                     "( 0  0   1 -15 )\n"
                                     "( 0  0   0   1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4   7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  0 -12  0 )
   // ( 0  0  1  1 )      ( 0  0   1  1 )
   // ( 0  0  0  1 )      ( 0  0   0  5 )
   {
      test_ = "Row-major sparse matrix subtraction assignment test 3";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1  -4   7 -2 )
   // ( 0  1  0  0 )  =>  ( 0   1 -12  0 )
   // ( 0  0  1  1 )      ( 0 -13   1  1 )
   // ( 0  0  0  1 )      ( 0   0   0  5 )
   {
      test_ = "Row-major sparse matrix subtraction assignment test 4";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 2UL );
      mat(0,1) = 12;
      mat(1,0) = 13;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1 -16  -8 -2 )
   // ( 0  1  0  0 )  =>  ( 0   1 -17  0 )
   // ( 0  0  1  1 )      ( 0   0   1  1 )
   // ( 0  0  0  1 )      ( 0   0   0  1 )
   {
      test_ = "Column-major sparse matrix subtraction assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 4UL );
      mat(0,0) = 12;
      mat(0,1) = 15;
      mat(1,1) = 17;
      mat.insert( 3UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 9UL );

      if( sm(0,0) != -16 || sm(0,1) !=  -8 ||
          sm(1,0) !=   1 || sm(1,1) != -17 ||
          sm(2,0) !=   0 || sm(2,1) !=   1 ||
          sm(3,0) !=   0 || sm(3,1) !=   0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( -16  -8 )\n(   1 -17 )\n(   0   1 )\n(   0   0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -16 || upper(0,2) !=  -8 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=   1 || upper(1,2) != -17 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=   0 || upper(2,2) !=   1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) !=   0 || upper(3,2) !=   0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -16  -8 -2 )\n"
                                     "( 0   1 -17  0 )\n"
                                     "( 0   0   1  1 )\n"
                                     "( 0   0   0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4   7  -2 )
   // ( 0  1  0  0 )  =>  ( 0  1 -14 -10 )
   // ( 0  0  1  1 )      ( 0  0   1 -15 )
   // ( 0  0  0  1 )      ( 0  0   0   1 )
   {
      test_ = "Column-major sparse matrix subtraction assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 4UL );
      mat(0,2) = 14;
      mat(0,3) = 10;
      mat(1,3) = 16;
      mat.insert( 1UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( upper,  4UL );
      checkColumns ( upper,  4UL );
      checkNonZeros( upper, 10UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != -14 || sm(0,3) != -10 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) !=   1 || sm(1,3) != -15 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1 -14 -10 )\n( 0  0   1 -15 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) !=   7 || upper(0,3) !=  -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != -14 || upper(1,3) != -10 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) !=   1 || upper(2,3) != -15 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) !=   0 || upper(3,3) !=   1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4   7  -2 )\n"
                                     "( 0  1 -14 -10 )\n"
                                     "( 0  0   1 -15 )\n"
                                     "( 0  0   0   1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4   7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  0 -12  0 )
   // ( 0  0  1  1 )      ( 0  0   1  1 )
   // ( 0  0  0  1 )      ( 0  0   0  5 )
   {
      test_ = "Column-major sparse matrix subtraction assignment test 3";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 12;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1  -4   7 -2 )
   // ( 0  1  0  0 )  =>  ( 0   1 -12  0 )
   // ( 0  0  1  1 )      ( 0 -13   1  1 )
   // ( 0  0  0  1 )      ( 0   0   0  5 )
   {
      test_ = "Column-major sparse matrix subtraction assignment test 4";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 2UL );
      mat(0,1) = 12;
      mat(1,0) = 13;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Schur product assignment to a submatrix of a UniUpperMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment to a submatrix of a UniUpperMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename UT >  // Type of the uniupper matrix
void SubmatrixTest::testSchurAssign()
{
   //=====================================================================================
   // Dense matrix Schur product assignment
   //=====================================================================================

   // ( 1 -4  7 -2 )      ( 1  8 21 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major dense matrix Schur product assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(0,0) = -2;
      mat(0,1) =  3;
      mat(1,0) =  1;
      mat(1,1) = 99;
      mat(2,1) =  1;
      mat(3,0) = 99;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 8UL );

      if( sm(0,0) != 8 || sm(0,1) != 21 ||
          sm(1,0) != 1 || sm(1,1) !=  0 ||
          sm(2,0) != 0 || sm(2,1) !=  1 ||
          sm(3,0) != 0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  8 21 )\n(  1  0 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 8 || upper(0,2) != 21 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) !=  0 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1  8 21 -2 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1 20 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major dense matrix Schur product assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) =  1;
      mat(0,2) = 99;
      mat(0,3) = 99;
      mat(1,0) = 99;
      mat(1,2) =  1;
      mat(1,3) = 20;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 8UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 0 || sm(0,3) !=  0 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) != 1 || sm(1,3) != 20 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1  0  0 )\n( 0  0  1 20 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 0 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 1 || upper(2,3) != 20 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) != 0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0  0  1 20 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  2  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Row-major dense matrix Schur product assignment test 3";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 99;
      mat(1,0) = 99;
      mat(1,1) =  2;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1  8 21 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major dense matrix Schur product assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(0,0) = -2;
      mat(0,1) =  3;
      mat(1,0) =  1;
      mat(1,1) = 99;
      mat(2,1) =  1;
      mat(3,0) = 99;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 8UL );

      if( sm(0,0) != 8 || sm(0,1) != 21 ||
          sm(1,0) != 1 || sm(1,1) !=  0 ||
          sm(2,0) != 0 || sm(2,1) !=  1 ||
          sm(3,0) != 0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  8 21 )\n(  1  0 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 8 || upper(0,2) != 21 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) !=  0 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1  8 21 -2 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1 20 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major dense matrix Schur product assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) =  1;
      mat(0,2) = 99;
      mat(0,3) = 99;
      mat(1,0) = 99;
      mat(1,2) =  1;
      mat(1,3) = 20;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 8UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 0 || sm(0,3) !=  0 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) != 1 || sm(1,3) != 20 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1  0  0 )\n( 0  0  1 20 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 0 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 1 || upper(2,3) != 20 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) != 0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0  0  1 20 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  2  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Column-major dense matrix Schur product assignment test 3";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL );
      mat(0,0) =  1;
      mat(0,1) = 99;
      mat(1,0) = 99;
      mat(1,1) =  2;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse matrix Schur product assignment
   //=====================================================================================

   // ( 1 -4  7 -2 )      ( 1  8 21 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major sparse matrix Schur product assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 7UL );
      mat(0,0) = -2;
      mat(0,1) =  3;
      mat(1,0) =  1;
      mat(1,1) = 99;
      mat(2,1) =  1;
      mat(3,0) = 99;
      mat.insert( 2UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 8UL );

      if( sm(0,0) != 8 || sm(0,1) != 21 ||
          sm(1,0) != 1 || sm(1,1) !=  0 ||
          sm(2,0) != 0 || sm(2,1) !=  1 ||
          sm(3,0) != 0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  8 21 )\n(  1  0 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 8 || upper(0,2) != 21 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) !=  0 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1  8 21 -2 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1 20 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Row-major sparse matrix Schur product assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 7UL );
      mat(0,1) =  1;
      mat(0,2) = 99;
      mat(0,3) = 99;
      mat(1,0) = 99;
      mat(1,2) =  1;
      mat(1,3) = 20;
      mat.insert( 1UL, 1UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 8UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 0 || sm(0,3) !=  0 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) != 1 || sm(1,3) != 20 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1  0  0 )\n( 0  0  1 20 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 0 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 1 || upper(2,3) != 20 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) != 0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0  0  1 20 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  2  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Row-major sparse matrix Schur product assignment test 3";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  1;
      mat(0,1) = 99;
      mat(1,0) = 99;
      mat(1,1) =  2;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1 -4  7 -2 )      ( 1  8 21 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major sparse matrix Schur product assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 7UL );
      mat(0,0) = -2;
      mat(0,1) =  3;
      mat(1,0) =  1;
      mat(1,1) = 99;
      mat(2,1) =  1;
      mat(3,0) = 99;
      mat.insert( 2UL, 0UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 8UL );

      if( sm(0,0) != 8 || sm(0,1) != 21 ||
          sm(1,0) != 1 || sm(1,1) !=  0 ||
          sm(2,0) != 0 || sm(2,1) !=  1 ||
          sm(3,0) != 0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  8 21 )\n(  1  0 )\n(  0  1 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != 8 || upper(0,2) != 21 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) != 1 || upper(1,2) !=  0 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) != 0 || upper(2,2) !=  1 || upper(2,3) !=  1 ||
          upper(3,0) != 0 || upper(3,1) != 0 || upper(3,2) !=  0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1  8 21 -2 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0  0  1  1 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  1  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1 20 )
   // ( 0  0  0  1 )      ( 0  0  0  1 )
   {
      test_ = "Column-major sparse matrix Schur product assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 7UL );
      mat(0,1) =  1;
      mat(0,2) = 99;
      mat(0,3) = 99;
      mat(1,0) = 99;
      mat(1,2) =  1;
      mat(1,3) = 20;
      mat.insert( 1UL, 1UL, 0 );

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( upper, 4UL );
      checkColumns ( upper, 4UL );
      checkNonZeros( upper, 8UL );

      if( sm(0,0) != 0 || sm(0,1) != 1 || sm(0,2) != 0 || sm(0,3) !=  0 ||
          sm(1,0) != 0 || sm(1,1) != 0 || sm(1,2) != 1 || sm(1,3) != 20 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0  1  0  0 )\n( 0  0  1 20 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( upper(0,0) != 1 || upper(0,1) != -4 || upper(0,2) != 7 || upper(0,3) != -2 ||
          upper(1,0) != 0 || upper(1,1) !=  1 || upper(1,2) != 0 || upper(1,3) !=  0 ||
          upper(2,0) != 0 || upper(2,1) !=  0 || upper(2,2) != 1 || upper(2,3) != 20 ||
          upper(3,0) != 0 || upper(3,1) !=  0 || upper(3,2) != 0 || upper(3,3) !=  1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << upper << "\n"
             << "   Expected result:\n( 1 -4  7 -2 )\n"
                                     "( 0  1  0  0 )\n"
                                     "( 0  0  1 20 )\n"
                                     "( 0  0  0  1 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1 -4  7 -2 )      ( 1 -4  7 -2 )
   // ( 0  1  0  0 )  =>  ( 0  2  0  0 )
   // ( 0  0  1  1 )      ( 0  0  1  1 )
   // ( 0  0  0  1 )      ( 0  0  0  5 )
   {
      test_ = "Column-major sparse matrix Schur product assignment test 3";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  1;
      mat(0,1) = 99;
      mat(1,0) = 99;
      mat(1,1) =  2;

      UT upper;
      init( upper );

      auto sm = submatrix( upper, 1UL, 1UL, 2UL, 2UL );

      try {
         sm %= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
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
/*!\brief Initializing the given upper unitriangular matrix.
//
// \return void
//
// This function is called before each test case to initialize the given upper unitriangular
// matrix.
*/
template< typename UT >
void SubmatrixTest::init( UT& upper )
{
   upper.resize( 4UL );
   upper(0,1) = -4;
   upper(0,2) =  7;
   upper(0,3) = -2;
   upper(1,2) =  0;
   upper(1,3) =  0;
   upper(2,3) =  1;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the assignment to a submatrix of a UniUpperMatrix.
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
/*!\brief Macro for the execution of the UniUpperMatrix submatrix test.
*/
#define RUN_UNIUPPERMATRIX_SUBMATRIX_TEST \
   blazetest::mathtest::adaptors::uniuppermatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace uniuppermatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
