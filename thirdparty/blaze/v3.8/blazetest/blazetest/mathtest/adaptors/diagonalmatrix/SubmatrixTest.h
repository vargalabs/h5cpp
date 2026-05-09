//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/diagonalmatrix/SubmatrixTest.h
//  \brief Header file for the DiagonalMatrix submatrix test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_DIAGONALMATRIX_SUBMATRIXTEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_DIAGONALMATRIX_SUBMATRIXTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/Submatrix.h>
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
/*!\brief Auxiliary class for assignment tests to a submatrix of a DiagonalMatrix.
//
// This class performs assignment tests to a submatrix of a DiagonalMatrix. It performs a series
// of both compile time as well as runtime tests.
*/
class SubmatrixTest
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
   template< typename DT > void testAssignment ();
   template< typename DT > void testAddAssign  ();
   template< typename DT > void testSubAssign  ();
   template< typename DT > void testSchurAssign();

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
/*!\brief Test of the assignment to a submatrix of a DiagonalMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to a submatrix of a DiagonalMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename DT >  // Type of the diagonal matrix
void SubmatrixTest::testAssignment()
{
   //=====================================================================================
   // Dense matrix assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major dense matrix assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,1) = 11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major dense matrix assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 18;
      mat(1,2) = 11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  1  5 )      ( 0  0  0  5 )
   {
      test_ = "Row-major dense matrix assignment test 3";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  1  5 )      ( 0  0  0  5 )
   {
      test_ = "Row-major dense matrix assignment test 4";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 18;
      mat(2,1) = 11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 18;
      mat(1,2) = 11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix assignment test 3";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix assignment test 4";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse matrix assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 4UL );
      mat(1,0) = 18;
      mat(2,1) = 11;
      mat.insert( 1UL, 1UL, 0 );
      mat.insert( 2UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 18;
      mat(1,2) = 11;
      mat.insert( 0UL, 2UL, 0 );
      mat.insert( 1UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix assignment test 3";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;
      mat.insert( 1UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix assignment test 4";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;
      mat.insert( 0UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 4UL );
      mat(1,0) = 18;
      mat(2,1) = 11;
      mat.insert( 1UL, 1UL, 0 );
      mat.insert( 2UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 18;
      mat(1,2) = 11;
      mat.insert( 0UL, 2UL, 0 );
      mat.insert( 1UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm = mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix assignment test 3";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(0,1) = 12;
      mat(1,1) = 14;
      mat.insert( 1UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix assignment test 4";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) = 11;
      mat(1,0) = 13;
      mat(1,1) = 14;
      mat.insert( 0UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm = mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the addition assignment to a submatrix of a DiagonalMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment to a submatrix of a DiagonalMatrix. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename DT >  // Type of the diagonal matrix
void SubmatrixTest::testAddAssign()
{
   //=====================================================================================
   // Dense matrix addition assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major dense matrix addition assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 16;
      mat(2,1) =  8;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major dense matrix addition assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 16;
      mat(1,2) =  8;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  1  5 )      ( 0  0  0  5 )
   {
      test_ = "Row-major dense matrix addition assignment test 3";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) =  9;
      mat(0,1) = 12;
      mat(1,1) = 11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  1  5 )      ( 0  0  0  5 )
   {
      test_ = "Row-major dense matrix addition assignment test 4";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) =  9;
      mat(1,0) = 13;
      mat(1,1) = 11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix addition assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = 16;
      mat(2,1) =  8;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix addition assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 16;
      mat(1,2) =  8;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix addition assignment test 3";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) =  9;
      mat(0,1) = 12;
      mat(1,1) = 11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix addition assignment test 4";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) =  9;
      mat(1,0) = 13;
      mat(1,1) = 11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse matrix addition assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix addition assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 4UL );
      mat(1,0) = 16;
      mat(2,1) =  8;
      mat.insert( 1UL, 1UL, 0 );
      mat.insert( 2UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix addition assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 16;
      mat(1,2) =  8;
      mat.insert( 0UL, 2UL, 0 );
      mat.insert( 1UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix addition assignment test 3";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  9;
      mat(0,1) = 12;
      mat(1,1) = 11;
      mat.insert( 1UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix addition assignment test 4";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  9;
      mat(1,0) = 13;
      mat(1,1) = 11;
      mat.insert( 0UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix addition assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 4UL );
      mat(1,0) = 16;
      mat(2,1) =  8;
      mat.insert( 1UL, 1UL, 0 );
      mat.insert( 2UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm += mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix addition assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 16;
      mat(1,2) =  8;
      mat.insert( 0UL, 2UL, 0 );
      mat.insert( 1UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm += mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix addition assignment test 3";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  9;
      mat(0,1) = 12;
      mat(1,1) = 11;
      mat.insert( 1UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix addition assignment test 4";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  9;
      mat(1,0) = 13;
      mat(1,1) = 11;
      mat.insert( 0UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm += mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the subtraction assignment to a submatrix of a DiagonalMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment to a submatrix of a DiagonalMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename DT >  // Type of the diagonal matrix
void SubmatrixTest::testSubAssign()
{
   //=====================================================================================
   // Dense matrix subtraction assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major dense matrix subtraction assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = -16;
      mat(2,1) =  -8;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major dense matrix subtraction assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = -16;
      mat(1,2) =  -8;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  1  5 )      ( 0  0  0  5 )
   {
      test_ = "Row-major dense matrix subtraction assignment test 3";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) =  -9;
      mat(0,1) = -12;
      mat(1,1) = -11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  1  5 )      ( 0  0  0  5 )
   {
      test_ = "Row-major dense matrix subtraction assignment test 4";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) =  -9;
      mat(1,0) = -13;
      mat(1,1) = -11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix subtraction assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(1,0) = -16;
      mat(2,1) =  -8;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix subtraction assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = -16;
      mat(1,2) =  -8;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix subtraction assignment test 3";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) =  -9;
      mat(0,1) = -12;
      mat(1,1) = -11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix subtraction assignment test 4";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 0 );
      mat(0,0) =  -9;
      mat(1,0) = -13;
      mat(1,1) = -11;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse matrix subtraction assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix subtraction assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 4UL );
      mat(1,0) = -16;
      mat(2,1) =  -8;
      mat.insert( 1UL, 1UL, 0 );
      mat.insert( 2UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix subtraction assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = -16;
      mat(1,2) =  -8;
      mat.insert( 0UL, 2UL, 0 );
      mat.insert( 1UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix subtraction assignment test 3";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  -9;
      mat(0,1) = -12;
      mat(1,1) = -11;
      mat.insert( 1UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix subtraction assignment test 4";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  -9;
      mat(1,0) = -13;
      mat(1,1) = -11;
      mat.insert( 0UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix subtraction assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 4UL );
      mat(1,0) = -16;
      mat(2,1) =  -8;
      mat.insert( 1UL, 1UL, 0 );
      mat.insert( 2UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm -= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
          sm(1,0) != 18 || sm(1,1) !=  0 ||
          sm(2,0) !=  0 || sm(2,1) != 11 ||
          sm(3,0) !=  0 || sm(3,1) !=  0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n(  0  0 )\n( 18  0 )\n(  0 11 )\n(  0  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 18  0  0 )
   // ( 0  0  3  0 )      ( 0  0 11  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix subtraction assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = -16;
      mat(1,2) =  -8;
      mat.insert( 0UL, 2UL, 0 );
      mat.insert( 1UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm -= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 11 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 18  0  0 )\n( 0  0 11  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 18 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 11 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 18  0  0 )\n"
                                     "( 0  0 11  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11 12  0 )
   // ( 0  0  3  0 )      ( 0  0 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix subtraction assignment test 3";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  -9;
      mat(0,1) = -12;
      mat(1,1) = -11;
      mat.insert( 1UL, 0UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 11  0  0 )
   // ( 0  0  3  0 )      ( 0 13 14  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix subtraction assignment test 4";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 2UL, 4UL );
      mat(0,0) =  -9;
      mat(1,0) = -13;
      mat(1,1) = -11;
      mat.insert( 0UL, 1UL, 0 );

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 1UL, 2UL, 2UL );

      try {
         sm -= mat;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid matrix succeeded\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Schur product assignment to a submatrix of a DiagonalMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment to a submatrix of a
// DiagonalMatrix. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename DT >  // Type of the diagonal matrix
void SubmatrixTest::testSchurAssign()
{
   //=====================================================================================
   // Dense matrix Schur product assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major dense matrix Schur product assignment test 1";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 0 );
      mat(0,1) = 9;
      mat(1,0) = 6;
      mat(2,1) = 4;
      mat(3,0) = 9;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

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

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major dense matrix Schur product assignment test 2";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 6;
      mat(0,3) = 9;
      mat(1,0) = 9;
      mat(1,2) = 4;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 12  0  0 )\n( 0  0 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix Schur product assignment test 1";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 0 );
      mat(0,1) = 9;
      mat(1,0) = 6;
      mat(2,1) = 4;
      mat(3,0) = 9;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

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

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major dense matrix Schur product assignment test 2";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 0 );
      mat(0,1) = 6;
      mat(0,3) = 9;
      mat(1,0) = 9;
      mat(1,2) = 4;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 12  0  0 )\n( 0  0 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse matrix Schur product assignment
   //=====================================================================================

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix Schur product assignment test 1";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 4UL, 2UL, 4UL );
      mat(0,1) = 9;
      mat(1,0) = 6;
      mat(2,1) = 4;
      mat(3,0) = 9;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

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

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Row-major sparse matrix Schur product assignment test 2";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 6;
      mat(0,3) = 9;
      mat(1,0) = 9;
      mat(1,2) = 4;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 12  0  0 )\n( 0  0 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix Schur product assignment test 1";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 4UL );
      mat(0,1) = 9;
      mat(1,0) = 6;
      mat(2,1) = 4;
      mat(3,0) = 9;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 0UL, 1UL, 4UL, 2UL );
      sm %= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

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

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( 1  0  0  0 )      ( 1  0  0  0 )
   // ( 0  2  0  0 )  =>  ( 0 12  0  0 )
   // ( 0  0  3  0 )      ( 0  0 12  0 )
   // ( 0  0  0  4 )      ( 0  0  0  4 )
   {
      test_ = "Column-major sparse matrix Schur product assignment test 2";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat( 2UL, 4UL, 4UL );
      mat(0,1) = 6;
      mat(0,3) = 9;
      mat(1,0) = 9;
      mat(1,2) = 4;

      DT diag;
      init( diag );

      auto sm = submatrix( diag, 1UL, 0UL, 2UL, 4UL );
      sm %= mat;

      checkRows    ( diag, 4UL );
      checkColumns ( diag, 4UL );
      checkNonZeros( diag, 4UL );

      if( sm(0,0) != 0 || sm(0,1) != 12 || sm(0,2) !=  0 || sm(0,3) != 0 ||
          sm(1,0) != 0 || sm(1,1) !=  0 || sm(1,2) != 12 || sm(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << sm << "\n"
             << "   Expected result:\n( 0 12  0  0 )\n( 0  0 12  0 )\n";
         throw std::runtime_error( oss.str() );
      }

      if( diag(0,0) != 1 || diag(0,1) !=  0 || diag(0,2) !=  0 || diag(0,3) != 0 ||
          diag(1,0) != 0 || diag(1,1) != 12 || diag(1,2) !=  0 || diag(1,3) != 0 ||
          diag(2,0) != 0 || diag(2,1) !=  0 || diag(2,2) != 12 || diag(2,3) != 0 ||
          diag(3,0) != 0 || diag(3,1) !=  0 || diag(3,2) !=  0 || diag(3,3) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to submatrix failed\n"
             << " Details:\n"
             << "   Result:\n" << diag << "\n"
             << "   Expected result:\n( 1  0  0  0 )\n"
                                     "( 0 12  0  0 )\n"
                                     "( 0  0 12  0 )\n"
                                     "( 0  0  0  4 )\n";
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
/*!\brief Initializing the given diagonal matrix.
//
// \return void
//
// This function is called before each test case to initialize the given diagonal matrix.
*/
template< typename DT >
void SubmatrixTest::init( DT& diag )
{
   diag.resize( 4UL );
   diag(0,0) = 1;
   diag(1,1) = 2;
   diag(2,2) = 3;
   diag(3,3) = 4;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the assignment to a submatrix of a DiagonalMatrix.
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
/*!\brief Macro for the execution of the DiagonalMatrix submatrix test.
*/
#define RUN_DIAGONALMATRIX_SUBMATRIX_TEST \
   blazetest::mathtest::adaptors::diagonalmatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace diagonalmatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
