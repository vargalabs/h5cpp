//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/hermitianmatrix/ColumnComplexTest.h
//  \brief Header file for the HermitianMatrix column complex test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_HERMITIANMATRIX_COLUMNCOMPLEXTEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_HERMITIANMATRIX_COLUMNCOMPLEXTEST_H_


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
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace hermitianmatrix {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for assignment tests to a single column of a HermitianMatrix.
//
// This class performs assignment tests to a single column of a HermitianMatrix with complex
// element type. It performs a series of both compile time as well as runtime tests.
*/
class ColumnComplexTest
{
 private:
   //**Type definitions****************************************************************************
   //! Complex element type.
   using cplx = blaze::complex<int>;

   //! Type of the dense Hermitian matrix.
   using DHT = blaze::HermitianMatrix< blaze::DynamicMatrix<cplx,blaze::rowMajor> >;

   //! Opposite dense Hermitian matrix type.
   using DOHT = DHT::OppositeType;

   //! Type of the sparse Hermitian matrix.
   using SHT = blaze::HermitianMatrix< blaze::CompressedMatrix<cplx,blaze::rowMajor> >;

   //! Opposite sparse Hermitian matrix type.
   using SOHT = SHT::OppositeType;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit ColumnComplexTest();
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
   template< typename HT > void testAssignment();
   template< typename HT > void testAddAssign ();
   template< typename HT > void testSubAssign ();
   template< typename HT > void testMultAssign();

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
   template< typename HT > void init( HT& herm );
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
/*!\brief Test of the assignment to columns of a HermitianMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to a single column of a HermitianMatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename HT >  // Type of the Hermitian matrix
void ColumnComplexTest::testAssignment()
{
   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 0) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Dense vector assignment test 1";

      blaze::DynamicVector<cplx,blaze::columnVector> vec( 3UL );
      vec[0] = cplx(2, 1);
      vec[1] = cplx(8, 0);
      vec[2] = cplx(4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );
      col1 = vec;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 9UL );

      if( col1[0] != cplx(2,1) || col1[1] != cplx(8,0) || col1[2] != cplx(4,-2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (2,1) (8,0) (4,-2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(2,-1) || herm(1,1) != cplx(8, 0) || herm(1,2) != cplx(4,2) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(4,-2) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 1) (7,3) )\n"
                                     "( (2,-1) (8, 0) (4,2) )\n"
                                     "( (7,-3) (4,-2) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 9) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Dense vector assignment test 2";

      blaze::DynamicVector<cplx,blaze::columnVector> vec( 3UL );
      vec[0] = cplx(2, 1);
      vec[1] = cplx(8, 9);
      vec[2] = cplx(4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 0) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Sparse vector assignment test 1";

      blaze::CompressedVector<cplx,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = cplx(2, 1);
      vec[1] = cplx(8, 0);
      vec[2] = cplx(4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );
      col1 = vec;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 9UL );

      if( col1[0] != cplx(2,1) || col1[1] != cplx(8,0) || col1[2] != cplx(4,-2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (2,1) (8,0) (4,-2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(2,-1) || herm(1,1) != cplx(8, 0) || herm(1,2) != cplx(4,2) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(4,-2) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 1) (7,3) )\n"
                                     "( (2,-1) (8, 0) (4,2) )\n"
                                     "( (7,-3) (4,-2) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 9) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Dense vector assignment test 2";

      blaze::CompressedVector<cplx,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = cplx(2, 1);
      vec[1] = cplx(8, 9);
      vec[2] = cplx(4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );

      try {
         col1 = vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the addition assignment to columns of a HermitianMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment to a single column of a HermitianMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename HT >  // Type of the Hermitian matrix
void ColumnComplexTest::testAddAssign()
{
   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 0) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Dense vector addition assignment test 1";

      blaze::DynamicVector<cplx,blaze::columnVector> vec( 3UL );
      vec[0] = cplx(6, 2);
      vec[1] = cplx(6, 0);
      vec[2] = cplx(4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );
      col1 += vec;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 9UL );

      if( col1[0] != cplx(2,1) || col1[1] != cplx(8,0) || col1[2] != cplx(4,-2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (2,1) (8,0) (4,-2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(2,-1) || herm(1,1) != cplx(8, 0) || herm(1,2) != cplx(4,2) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(4,-2) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 1) (7,3) )\n"
                                     "( (2,-1) (8, 0) (4,2) )\n"
                                     "( (7,-3) (4,-2) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 9) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Dense vector addition assignment test 2";

      blaze::DynamicVector<cplx,blaze::columnVector> vec( 3UL );
      vec[0] = cplx(6, 2);
      vec[1] = cplx(6, 9);
      vec[2] = cplx(4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );

      try {
         col1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 0) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Sparse vector addition assignment test 1";

      blaze::CompressedVector<cplx,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = cplx(6, 2);
      vec[1] = cplx(6, 0);
      vec[2] = cplx(4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );
      col1 += vec;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 9UL );

      if( col1[0] != cplx(2,1) || col1[1] != cplx(8,0) || col1[2] != cplx(4,-2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (2,1) (8,0) (4,-2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(2,-1) || herm(1,1) != cplx(8, 0) || herm(1,2) != cplx(4,2) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(4,-2) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 1) (7,3) )\n"
                                     "( (2,-1) (8, 0) (4,2) )\n"
                                     "( (7,-3) (4,-2) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 9) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Sparse vector addition assignment test 2";

      blaze::CompressedVector<cplx,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = cplx(6, 2);
      vec[1] = cplx(6, 9);
      vec[2] = cplx(4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );

      try {
         col1 += vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the subtraction assignment to columns of a HermitianMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment to a single column of a
// HermitianMatrix. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename HT >  // Type of the Hermitian matrix
void ColumnComplexTest::testSubAssign()
{
   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 0) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Dense vector subtraction assignment test 1";

      blaze::DynamicVector<cplx,blaze::columnVector> vec( 3UL );
      vec[0] = cplx(-6,-2);
      vec[1] = cplx(-6, 0);
      vec[2] = cplx(-4, 2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );
      col1 -= vec;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 9UL );

      if( col1[0] != cplx(2,1) || col1[1] != cplx(8,0) || col1[2] != cplx(4,-2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (2,1) (8,0) (4,-2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(2,-1) || herm(1,1) != cplx(8, 0) || herm(1,2) != cplx(4,2) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(4,-2) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 1) (7,3) )\n"
                                     "( (2,-1) (8, 0) (4,2) )\n"
                                     "( (7,-3) (4,-2) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 9) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Dense vector subtraction assignment test 2";

      blaze::DynamicVector<cplx,blaze::columnVector> vec( 3UL );
      vec[0] = cplx(-6,-2);
      vec[1] = cplx(-6,-9);
      vec[2] = cplx(-4, 2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );

      try {
         col1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector subtraction assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 0) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Sparse vector subtraction assignment test 1";

      blaze::CompressedVector<cplx,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = cplx(-6,-2);
      vec[1] = cplx(-6, 0);
      vec[2] = cplx(-4, 2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );
      col1 -= vec;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 9UL );

      if( col1[0] != cplx(2,1) || col1[1] != cplx(8,0) || col1[2] != cplx(4,-2) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (2,1) (8,0) (4,-2) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(2, 1) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(2,-1) || herm(1,1) != cplx(8, 0) || herm(1,2) != cplx(4,2) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(4,-2) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (2, 1) (7,3) )\n"
                                     "( (2,-1) (8, 0) (4,2) )\n"
                                     "( (7,-3) (4,-2) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (2, 1) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (2,-1) (8, 9) (4,2) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (4,-2) (3,0) )
   {
      test_ = "Sparse vector subtraction assignment test 2";

      blaze::CompressedVector<cplx,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = cplx(-6,-2);
      vec[1] = cplx(-6,-9);
      vec[2] = cplx(-4, 2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );

      try {
         col1 -= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the multiplication assignment to columns of a HermitianMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment to a single column of a
// HermitianMatrix. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename HT >  // Type of the Hermitian matrix
void ColumnComplexTest::testMultAssign()
{
   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (8,2) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (8,-2) (6,0) (0,0) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (0,0) (3,0) )
   {
      test_ = "Dense vector multiplication assignment test 1";

      blaze::DynamicVector<cplx,blaze::columnVector> vec( 3UL );
      vec[0] = cplx(-2, 0);
      vec[1] = cplx( 3, 0);
      vec[2] = cplx( 4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );
      col1 *= vec;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( col1[0] != cplx(8,2) || col1[1] != cplx(6,0) || col1[2] != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (8,2) (6,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(8,2) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(8,-2) || herm(1,1) != cplx(6,0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (8,2) (7,3) )\n"
                                     "( (8,-2) (6,0) (0,0) )\n"
                                     "( (7,-3) (0,0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (8,2) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (8,-2) (6,4) (0,0) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (0,0) (3,0) )
   {
      test_ = "Dense vector multiplication assignment test 2";

      blaze::DynamicVector<cplx,blaze::columnVector> vec( 3UL );
      vec[0] = cplx(-2, 0);
      vec[1] = cplx( 3, 2);
      vec[2] = cplx( 4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );

      try {
         col1 *= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::invalid_argument& ) {}
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (8,2) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (8,-2) (6,0) (0,0) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (0,0) (3,0) )
   {
      test_ = "Sparse vector multiplication assignment test 1";

      blaze::CompressedVector<cplx,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = cplx(-2, 0);
      vec[1] = cplx( 3, 0);
      vec[2] = cplx( 4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );
      col1 *= vec;

      checkRows    ( herm, 3UL );
      checkColumns ( herm, 3UL );
      checkNonZeros( herm, 7UL );

      if( col1[0] != cplx(8,2) || col1[1] != cplx(6,0) || col1[2] != cplx(0,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( (8,2) (6,0) (0,0) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( herm(0,0) != cplx(1, 0) || herm(0,1) != cplx(8,2) || herm(0,2) != cplx(7,3) ||
          herm(1,0) != cplx(8,-2) || herm(1,1) != cplx(6,0) || herm(1,2) != cplx(0,0) ||
          herm(2,0) != cplx(7,-3) || herm(2,1) != cplx(0,0) || herm(2,2) != cplx(3,0) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment to column failed\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n"
             << "   Expected result:\n( (1, 0) (8,2) (7,3) )\n"
                                     "( (8,-2) (6,0) (0,0) )\n"
                                     "( (7,-3) (0,0) (3,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // ( ( 1, 0) (-4,-1) (7,3) )      ( (1, 0) (8,2) (7,3) )
   // ( (-4, 1) ( 2, 0) (0,0) )  =>  ( (8,-2) (6,4) (0,0) )
   // ( ( 7,-3) ( 0, 0) (3,0) )      ( (7,-3) (0,0) (3,0) )
   {
      test_ = "Sparse vector multiplication assignment test 2";

      blaze::CompressedVector<cplx,blaze::columnVector> vec( 3UL, 3UL );
      vec[0] = cplx(-2, 0);
      vec[1] = cplx( 3, 2);
      vec[2] = cplx( 4,-2);

      HT herm;
      init( herm );

      auto col1 = column( herm, 1UL );

      try {
         col1 *= vec;

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment of invalid vector succeeded\n"
             << " Details:\n"
             << "   Result:\n" << herm << "\n";
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
void ColumnComplexTest::checkRows( const Type& matrix, size_t expectedRows ) const
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
void ColumnComplexTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
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
void ColumnComplexTest::checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const
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
/*!\brief Initializing the given Hermitian matrix.
//
// \return void
//
// This function is called before each test case to initialize the given Hermitian matrix.
*/
template< typename HT >
void ColumnComplexTest::init( HT& herm )
{
   herm.resize( 3UL );
   herm(0,0) = cplx( 1, 0);
   herm(0,1) = cplx(-4,-1);
   herm(0,2) = cplx( 7, 3);
   herm(1,1) = cplx( 2, 0);
   herm(2,2) = cplx( 3, 0);
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the assignment to a single column of a HermitianMatrix.
//
// \return void
*/
void runTest()
{
   ColumnComplexTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the HermitianMatrix column complex test.
*/
#define RUN_HERMITIANMATRIX_COLUMNCOMPLEX_TEST \
   blazetest::mathtest::adaptors::hermitianmatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace hermitianmatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
