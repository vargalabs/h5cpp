//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/symmetricmatrix/SubmatrixScalarTest.h
//  \brief Header file for the SymmetricMatrix submatrix scalar test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_SYMMETRICMATRIX_SUBMATRIXSCALARTEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_SYMMETRICMATRIX_SUBMATRIXSCALARTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Submatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace symmetricmatrix {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for assignment tests to a submatrix of a scalar SymmetricMatrix.
//
// This class performs assignment tests to a submatrix of a SymmetricMatrix with scalar element
// type. It performs a series of both compile time as well as runtime tests.
*/
class SubmatrixScalarTest
{
 private:
   //**Type definitions****************************************************************************
   //! Type of the dense scalar symmetric matrix.
   using DST = blaze::SymmetricMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> >;

   //! Opposite dense scalar symmetric matrix type.
   using DOST = DST::OppositeType;

   //! Type of the sparse scalar symmetric matrix.
   using SST = blaze::SymmetricMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> >;

   //! Opposite sparse scalar symmetric matrix type.
   using SOST = SST::OppositeType;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit SubmatrixScalarTest();
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
   template< typename ST > void testAssignment ();
   template< typename ST > void testAddAssign  ();
   template< typename ST > void testSubAssign  ();
   template< typename ST > void testSchurAssign();

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
   template< typename ST > void init( ST& sym );
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
/*!\brief Test of the assignment to a submatrix of a SymmetricMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to a submatrix of a SymmetricMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename ST >  // Type of the symmetric matrix
void SubmatrixScalarTest::testAssignment()
{
   //=====================================================================================
   // Dense matrix assignment
   //=====================================================================================

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 18 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix assignment test 1";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 15;
         mat(1,0) = 18;
         mat(1,1) = 17;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 18 || sm(1,1) != 17 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 18 17 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(1,0) = 18;
         mat(1,1) = 17;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(3,0) = 15;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 ||
             sm(1,0) != 18 || sm(1,1) != 17 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 )\n"
                                        "( 18 17 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 14 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix assignment test 2";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 15;
         mat(1,0) = 13;
         mat(1,1) = 14;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(3,0) = 15;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 11 19 )
   {
      test_ = "Dense matrix assignment test 3";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 11;
         mat(1,0) = 13;
         mat(1,1) = 14;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(3,0) = 11;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 11 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 11 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix assignment test 4";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 11;
         mat(1,0) = 13;
         mat(1,1) = 14;
         mat(1,2) = 11;
         mat(1,3) = 19;
         mat(2,0) = 19;
         mat(2,1) = 11;
         mat(2,2) = 12;
         mat(2,3) = 14;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ||
             sm(2,0) != 19 || sm(2,1) != 11 || sm(2,2) != 12 || sm(2,3) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n"
                                        "( 13 14 11 19 )\n"
                                        "( 19 11 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(0,2) = 19;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(1,2) = 11;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(2,2) = 12;
         mat(3,0) = 11;
         mat(3,1) = 19;
         mat(3,2) = 14;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 || sm(0,2) != 19 ||
             sm(1,0) != 18 || sm(1,1) != 14 || sm(1,2) != 11 ||
             sm(2,0) != 14 || sm(2,1) != 11 || sm(2,2) != 12 ||
             sm(3,0) != 11 || sm(3,1) != 19 || sm(3,2) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 19 )\n"
                                        "( 18 14 11 )\n"
                                        "( 14 11 12 )\n"
                                        "( 11 19 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 22 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix assignment test 5";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 15;
         mat(1,0) = 22;
         mat(1,1) = 17;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(1,0) = 22;
         mat(1,1) = 17;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(3,0) = 15;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 22 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix assignment test 6";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 15;
         mat(1,0) = 13;
         mat(1,1) = 22;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(2,0) = 22;
         mat(2,1) = 11;
         mat(3,0) = 15;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 22 19 )
   {
      test_ = "Dense matrix assignment test 7";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 11;
         mat(1,0) = 13;
         mat(1,1) = 14;
         mat(1,2) = 22;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(3,0) = 22;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 22 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix assignment test 8";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 11;
         mat(1,0) = 13;
         mat(1,1) = 14;
         mat(1,2) = 22;
         mat(1,3) = 19;
         mat(2,0) = 19;
         mat(2,1) = 11;
         mat(2,2) = 12;
         mat(2,3) = 14;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(0,2) = 19;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(1,2) = 11;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(2,2) = 12;
         mat(3,0) = 22;
         mat(3,1) = 19;
         mat(3,2) = 14;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Sparse matrix assignment
   //=====================================================================================

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 18 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix assignment test 1";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 15;
         mat(1,0) = 18;
         mat(1,1) = 17;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 18 || sm(1,1) != 17 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 18 17 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(1,0) = 18;
         mat(1,1) = 17;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(3,0) = 15;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 ||
             sm(1,0) != 18 || sm(1,1) != 17 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 )\n"
                                        "( 18 17 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 14 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix assignment test 2";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 15;
         mat(1,0) = 13;
         mat(1,1) = 14;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(3,0) = 15;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 11 19 )
   {
      test_ = "Sparse matrix assignment test 3";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 11;
         mat(1,0) = 13;
         mat(1,1) = 14;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(3,0) = 11;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 11 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 11 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix assignment test 4";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 11;
         mat(1,0) = 13;
         mat(1,1) = 14;
         mat(1,2) = 11;
         mat(1,3) = 19;
         mat(2,0) = 19;
         mat(2,1) = 11;
         mat(2,2) = 12;
         mat(2,3) = 14;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ||
             sm(2,0) != 19 || sm(2,1) != 11 || sm(2,2) != 12 || sm(2,3) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n"
                                        "( 13 14 11 19 )\n"
                                        "( 19 11 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(0,2) = 19;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(1,2) = 11;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(2,2) = 12;
         mat(3,0) = 11;
         mat(3,1) = 19;
         mat(3,2) = 14;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );
         sm = mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 || sm(0,2) != 19 ||
             sm(1,0) != 18 || sm(1,1) != 14 || sm(1,2) != 11 ||
             sm(2,0) != 14 || sm(2,1) != 11 || sm(2,2) != 12 ||
             sm(3,0) != 11 || sm(3,1) != 19 || sm(3,2) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 19 )\n"
                                        "( 18 14 11 )\n"
                                        "( 14 11 12 )\n"
                                        "( 11 19 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 22 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix assignment test 5";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 15;
         mat(1,0) = 22;
         mat(1,1) = 17;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(1,0) = 22;
         mat(1,1) = 17;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(3,0) = 15;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 22 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix assignment test 6";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 15;
         mat(1,0) = 13;
         mat(1,1) = 22;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(2,0) = 22;
         mat(2,1) = 11;
         mat(3,0) = 15;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 22 19 )
   {
      test_ = "Sparse matrix assignment test 7";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 11;
         mat(1,0) = 13;
         mat(1,1) = 14;
         mat(1,2) = 22;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(3,0) = 22;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 22 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix assignment test 8";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = 12;
         mat(0,1) = 18;
         mat(0,2) = 14;
         mat(0,3) = 11;
         mat(1,0) = 13;
         mat(1,1) = 14;
         mat(1,2) = 22;
         mat(1,3) = 19;
         mat(2,0) = 19;
         mat(2,1) = 11;
         mat(2,2) = 12;
         mat(2,3) = 14;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(0,2) = 19;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(1,2) = 11;
         mat(2,0) = 14;
         mat(2,1) = 11;
         mat(2,2) = 12;
         mat(3,0) = 22;
         mat(3,1) = 19;
         mat(3,2) = 14;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the addition assignment to a submatrix of a SymmetricMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment to a submatrix of a SymmetricMatrix. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename ST >  // Type of the symmetric matrix
void SubmatrixScalarTest::testAddAssign()
{
   //=====================================================================================
   // Dense matrix addition assignment
   //=====================================================================================

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 18 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix addition assignment test 1";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 11;
         mat(0,1) = 22;
         mat(0,2) =  7;
         mat(0,3) = 17;
         mat(1,0) = 22;
         mat(1,1) = 15;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 18 || sm(1,1) != 17 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 18 17 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 11;
         mat(0,1) = 22;
         mat(1,0) = 22;
         mat(1,1) = 15;
         mat(2,0) =  7;
         mat(2,1) = 11;
         mat(3,0) = 17;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 ||
             sm(1,0) != 18 || sm(1,1) != 17 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 )\n"
                                        "( 18 17 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 14 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix addition assignment test 2";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 15;
         mat(0,2) = 13;
         mat(0,3) = 15;
         mat(1,0) = 13;
         mat(1,1) = 13;
         mat(1,2) =  6;
         mat(1,3) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 15;
         mat(1,1) = 13;
         mat(2,0) = 13;
         mat(2,1) =  6;
         mat(3,0) = 15;
         mat(3,1) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 11 19 )
   {
      test_ = "Dense matrix addition assignment test 3";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 11;
         mat(0,2) = 13;
         mat(0,3) = 15;
         mat(1,0) = 15;
         mat(1,1) = 14;
         mat(1,2) = 15;
         mat(1,3) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 12;
         mat(0,1) = 15;
         mat(1,0) = 11;
         mat(1,1) = 14;
         mat(2,0) = 13;
         mat(2,1) = 15;
         mat(3,0) = 15;
         mat(3,1) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 11 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 11 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix addition assignment test 4";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) =  5;
         mat(0,1) = 18;
         mat(0,2) = 11;
         mat(0,3) = 10;
         mat(1,0) = 15;
         mat(1,1) = 14;
         mat(1,2) = 10;
         mat(1,3) = 14;
         mat(2,0) = 14;
         mat(2,1) = 12;
         mat(2,2) = 12;
         mat(2,3) =  7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ||
             sm(2,0) != 19 || sm(2,1) != 11 || sm(2,2) != 12 || sm(2,3) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n"
                                        "( 13 14 11 19 )\n"
                                        "( 19 11 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) =  5;
         mat(0,1) = 15;
         mat(0,2) = 14;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(1,2) = 12;
         mat(2,0) = 11;
         mat(2,1) = 10;
         mat(2,2) = 12;
         mat(3,0) = 10;
         mat(3,1) = 14;
         mat(3,2) =  7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 || sm(0,2) != 19 ||
             sm(1,0) != 18 || sm(1,1) != 14 || sm(1,2) != 11 ||
             sm(2,0) != 14 || sm(2,1) != 11 || sm(2,2) != 12 ||
             sm(3,0) != 11 || sm(3,1) != 19 || sm(3,2) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 19 )\n"
                                        "( 18 14 11 )\n"
                                        "( 14 11 12 )\n"
                                        "( 11 19 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 22 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix addition assignment test 5";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 11;
         mat(0,1) = 22;
         mat(0,2) =  7;
         mat(0,3) = 17;
         mat(1,0) = 26;
         mat(1,1) = 15;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 11;
         mat(0,1) = 22;
         mat(1,0) = 26;
         mat(1,1) = 15;
         mat(2,0) =  7;
         mat(2,1) = 11;
         mat(3,0) = 17;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 22 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix addition assignment test 6";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 15;
         mat(0,2) = 13;
         mat(0,3) = 15;
         mat(1,0) = 13;
         mat(1,1) = 21;
         mat(1,2) =  6;
         mat(1,3) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 15;
         mat(1,1) = 13;
         mat(2,0) = 21;
         mat(2,1) =  6;
         mat(3,0) = 15;
         mat(3,1) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 22 19 )
   {
      test_ = "Dense matrix addition assignment test 7";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 12;
         mat(0,1) = 11;
         mat(0,2) = 13;
         mat(0,3) = 15;
         mat(1,0) = 15;
         mat(1,1) = 14;
         mat(1,2) = 26;
         mat(1,3) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 12;
         mat(0,1) = 15;
         mat(1,0) = 11;
         mat(1,1) = 14;
         mat(2,0) = 13;
         mat(2,1) = 15;
         mat(3,0) = 26;
         mat(3,1) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 22 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix addition assignment test 8";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) =  5;
         mat(0,1) = 18;
         mat(0,2) = 11;
         mat(0,3) = 10;
         mat(1,0) = 15;
         mat(1,1) = 14;
         mat(1,2) = 21;
         mat(1,3) = 14;
         mat(2,0) = 14;
         mat(2,1) = 12;
         mat(2,2) = 12;
         mat(2,3) =  7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) =  5;
         mat(0,1) = 15;
         mat(0,2) = 14;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(1,2) = 12;
         mat(2,0) = 11;
         mat(2,1) = 10;
         mat(2,2) = 12;
         mat(3,0) = 21;
         mat(3,1) = 14;
         mat(3,2) =  7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Sparse matrix addition assignment
   //=====================================================================================

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 18 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix addition assignment test 1";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 11;
         mat(0,1) = 22;
         mat(0,2) =  7;
         mat(0,3) = 17;
         mat(1,0) = 22;
         mat(1,1) = 15;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 18 || sm(1,1) != 17 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 18 17 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 11;
         mat(0,1) = 22;
         mat(1,0) = 22;
         mat(1,1) = 15;
         mat(2,0) =  7;
         mat(2,1) = 11;
         mat(3,0) = 17;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 ||
             sm(1,0) != 18 || sm(1,1) != 17 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 )\n"
                                        "( 18 17 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 14 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix addition assignment test 2";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 15;
         mat(0,2) = 13;
         mat(0,3) = 15;
         mat(1,0) = 13;
         mat(1,1) = 13;
         mat(1,2) =  6;
         mat(1,3) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 15;
         mat(1,1) = 13;
         mat(2,0) = 13;
         mat(2,1) =  6;
         mat(3,0) = 15;
         mat(3,1) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 11 19 )
   {
      test_ = "Sparse matrix addition assignment test 3";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 11;
         mat(0,2) = 13;
         mat(0,3) = 15;
         mat(1,0) = 15;
         mat(1,1) = 14;
         mat(1,2) = 15;
         mat(1,3) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 15;
         mat(1,0) = 11;
         mat(1,1) = 14;
         mat(2,0) = 13;
         mat(2,1) = 15;
         mat(3,0) = 15;
         mat(3,1) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 11 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 11 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix addition assignment test 4";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) =  5;
         mat(0,1) = 18;
         mat(0,2) = 11;
         mat(0,3) = 10;
         mat(1,0) = 15;
         mat(1,1) = 14;
         mat(1,2) = 10;
         mat(1,3) = 14;
         mat(2,0) = 14;
         mat(2,1) = 12;
         mat(2,2) = 12;
         mat(2,3) =  7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ||
             sm(2,0) != 19 || sm(2,1) != 11 || sm(2,2) != 12 || sm(2,3) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n"
                                        "( 13 14 11 19 )\n"
                                        "( 19 11 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) =  5;
         mat(0,1) = 15;
         mat(0,2) = 14;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(1,2) = 12;
         mat(2,0) = 11;
         mat(2,1) = 10;
         mat(2,2) = 12;
         mat(3,0) = 10;
         mat(3,1) = 14;
         mat(3,2) =  7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );
         sm += mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 || sm(0,2) != 19 ||
             sm(1,0) != 18 || sm(1,1) != 14 || sm(1,2) != 11 ||
             sm(2,0) != 14 || sm(2,1) != 11 || sm(2,2) != 12 ||
             sm(3,0) != 11 || sm(3,1) != 19 || sm(3,2) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 19 )\n"
                                        "( 18 14 11 )\n"
                                        "( 14 11 12 )\n"
                                        "( 11 19 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 22 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix addition assignment test 5";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 11;
         mat(0,1) = 22;
         mat(0,2) =  7;
         mat(0,3) = 17;
         mat(1,0) = 26;
         mat(1,1) = 15;
         mat(1,2) = 11;
         mat(1,3) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 11;
         mat(0,1) = 22;
         mat(1,0) = 26;
         mat(1,1) = 15;
         mat(2,0) =  7;
         mat(2,1) = 11;
         mat(3,0) = 17;
         mat(3,1) = 19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 22 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix addition assignment test 6";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 15;
         mat(0,2) = 13;
         mat(0,3) = 15;
         mat(1,0) = 13;
         mat(1,1) = 21;
         mat(1,2) =  6;
         mat(1,3) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 13;
         mat(1,0) = 15;
         mat(1,1) = 13;
         mat(2,0) = 21;
         mat(2,1) =  6;
         mat(3,0) = 15;
         mat(3,1) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 22 19 )
   {
      test_ = "Sparse matrix addition assignment test 7";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 11;
         mat(0,2) = 13;
         mat(0,3) = 15;
         mat(1,0) = 15;
         mat(1,1) = 14;
         mat(1,2) = 26;
         mat(1,3) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 12;
         mat(0,1) = 15;
         mat(1,0) = 11;
         mat(1,1) = 14;
         mat(2,0) = 13;
         mat(2,1) = 15;
         mat(3,0) = 26;
         mat(3,1) = 12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 22 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix addition assignment test 8";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) =  5;
         mat(0,1) = 18;
         mat(0,2) = 11;
         mat(0,3) = 10;
         mat(1,0) = 15;
         mat(1,1) = 14;
         mat(1,2) = 21;
         mat(1,3) = 14;
         mat(2,0) = 14;
         mat(2,1) = 12;
         mat(2,2) = 12;
         mat(2,3) =  7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) =  5;
         mat(0,1) = 15;
         mat(0,2) = 14;
         mat(1,0) = 18;
         mat(1,1) = 14;
         mat(1,2) = 12;
         mat(2,0) = 11;
         mat(2,1) = 10;
         mat(2,2) = 12;
         mat(3,0) = 21;
         mat(3,1) = 14;
         mat(3,2) =  7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the subtraction assignment to a submatrix of a SymmetricMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment to a submatrix of a SymmetricMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename ST >  // Type of the symmetric matrix
void SubmatrixScalarTest::testSubAssign()
{
   //=====================================================================================
   // Dense matrix subtraction assignment
   //=====================================================================================

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 18 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix subtraction assignment test 1";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = -11;
         mat(0,1) = -22;
         mat(0,2) =  -7;
         mat(0,3) = -17;
         mat(1,0) = -22;
         mat(1,1) = -15;
         mat(1,2) = -11;
         mat(1,3) = -19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 18 || sm(1,1) != 17 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 18 17 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = -11;
         mat(0,1) = -22;
         mat(1,0) = -22;
         mat(1,1) = -15;
         mat(2,0) =  -7;
         mat(2,1) = -11;
         mat(3,0) = -17;
         mat(3,1) = -19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 ||
             sm(1,0) != 18 || sm(1,1) != 17 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 )\n"
                                        "( 18 17 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 14 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix subtraction assignment test 2";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = -12;
         mat(0,1) = -15;
         mat(0,2) = -13;
         mat(0,3) = -15;
         mat(1,0) = -13;
         mat(1,1) = -13;
         mat(1,2) =  -6;
         mat(1,3) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = -12;
         mat(0,1) = -13;
         mat(1,0) = -15;
         mat(1,1) = -13;
         mat(2,0) = -13;
         mat(2,1) =  -6;
         mat(3,0) = -15;
         mat(3,1) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 11 19 )
   {
      test_ = "Dense matrix subtraction assignment test 3";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = -12;
         mat(0,1) = -11;
         mat(0,2) = -13;
         mat(0,3) = -15;
         mat(1,0) = -15;
         mat(1,1) = -14;
         mat(1,2) = -15;
         mat(1,3) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = -12;
         mat(0,1) = -15;
         mat(1,0) = -11;
         mat(1,1) = -14;
         mat(2,0) = -13;
         mat(2,1) = -15;
         mat(3,0) = -15;
         mat(3,1) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 11 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 11 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix subtraction assignment test 4";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) =  -5;
         mat(0,1) = -18;
         mat(0,2) = -11;
         mat(0,3) = -10;
         mat(1,0) = -15;
         mat(1,1) = -14;
         mat(1,2) = -10;
         mat(1,3) = -14;
         mat(2,0) = -14;
         mat(2,1) = -12;
         mat(2,2) = -12;
         mat(2,3) =  -7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ||
             sm(2,0) != 19 || sm(2,1) != 11 || sm(2,2) != 12 || sm(2,3) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n"
                                        "( 13 14 11 19 )\n"
                                        "( 19 11 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) =  -5;
         mat(0,1) = -15;
         mat(0,2) = -14;
         mat(1,0) = -18;
         mat(1,1) = -14;
         mat(1,2) = -12;
         mat(2,0) = -11;
         mat(2,1) = -10;
         mat(2,2) = -12;
         mat(3,0) = -10;
         mat(3,1) = -14;
         mat(3,2) =  -7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 || sm(0,2) != 19 ||
             sm(1,0) != 18 || sm(1,1) != 14 || sm(1,2) != 11 ||
             sm(2,0) != 14 || sm(2,1) != 11 || sm(2,2) != 12 ||
             sm(3,0) != 11 || sm(3,1) != 19 || sm(3,2) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 19 )\n"
                                        "( 18 14 11 )\n"
                                        "( 14 11 12 )\n"
                                        "( 11 19 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 22 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix subtraction assignment test 5";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = -11;
         mat(0,1) = -22;
         mat(0,2) =  -7;
         mat(0,3) = -17;
         mat(1,0) = -26;
         mat(1,1) = -15;
         mat(1,2) = -11;
         mat(1,3) = -19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = -11;
         mat(0,1) = -22;
         mat(1,0) = -26;
         mat(1,1) = -15;
         mat(2,0) =  -7;
         mat(2,1) = -11;
         mat(3,0) = -17;
         mat(3,1) = -19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 22 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix subtraction assignment test 6";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = -12;
         mat(0,1) = -15;
         mat(0,2) = -13;
         mat(0,3) = -15;
         mat(1,0) = -13;
         mat(1,1) = -21;
         mat(1,2) =  -6;
         mat(1,3) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = -12;
         mat(0,1) = -13;
         mat(1,0) = -15;
         mat(1,1) = -13;
         mat(2,0) = -21;
         mat(2,1) =  -6;
         mat(3,0) = -15;
         mat(3,1) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 22 19 )
   {
      test_ = "Dense matrix subtraction assignment test 7";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = -12;
         mat(0,1) = -11;
         mat(0,2) = -13;
         mat(0,3) = -15;
         mat(1,0) = -15;
         mat(1,1) = -14;
         mat(1,2) = -26;
         mat(1,3) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = -12;
         mat(0,1) = -15;
         mat(1,0) = -11;
         mat(1,1) = -14;
         mat(2,0) = -13;
         mat(2,1) = -15;
         mat(3,0) = -26;
         mat(3,1) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 22 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix subtraction assignment test 8";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) =  -5;
         mat(0,1) = -18;
         mat(0,2) = -11;
         mat(0,3) = -10;
         mat(1,0) = -15;
         mat(1,1) = -14;
         mat(1,2) = -21;
         mat(1,3) = -14;
         mat(2,0) = -14;
         mat(2,1) = -12;
         mat(2,2) = -12;
         mat(2,3) =  -7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) =  -5;
         mat(0,1) = -15;
         mat(0,2) = -14;
         mat(1,0) = -18;
         mat(1,1) = -14;
         mat(1,2) = -12;
         mat(2,0) = -11;
         mat(2,1) = -10;
         mat(2,2) = -12;
         mat(3,0) = -21;
         mat(3,1) = -14;
         mat(3,2) =  -7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Sparse matrix subtraction assignment
   //=====================================================================================

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 18 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix subtraction assignment test 1";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = -11;
         mat(0,1) = -22;
         mat(0,2) =  -7;
         mat(0,3) = -17;
         mat(1,0) = -22;
         mat(1,1) = -15;
         mat(1,2) = -11;
         mat(1,3) = -19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 18 || sm(1,1) != 17 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 18 17 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = -11;
         mat(0,1) = -22;
         mat(1,0) = -22;
         mat(1,1) = -15;
         mat(2,0) =  -7;
         mat(2,1) = -11;
         mat(3,0) = -17;
         mat(3,1) = -19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 ||
             sm(1,0) != 18 || sm(1,1) != 17 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 )\n"
                                        "( 18 17 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 12 || sym(0,1) != 18 || sym(0,2) != 14 || sym(0,3) != 15 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 18 || sym(1,1) != 17 || sym(1,2) != 11 || sym(1,3) != 19 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) != 11 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 15 || sym(3,1) != 19 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 12 18 14 15  5  0 )\n"
                                        "( 18 17 11 19 -1  8 )\n"
                                        "( 14 11  3  1  0 -2 )\n"
                                        "( 15 19  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 14 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix subtraction assignment test 2";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = -12;
         mat(0,1) = -15;
         mat(0,2) = -13;
         mat(0,3) = -15;
         mat(1,0) = -13;
         mat(1,1) = -13;
         mat(1,2) =  -6;
         mat(1,3) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 15 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 15 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = -12;
         mat(0,1) = -13;
         mat(1,0) = -15;
         mat(1,1) = -13;
         mat(2,0) = -13;
         mat(2,1) =  -6;
         mat(3,0) = -15;
         mat(3,1) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 15 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 15 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 12 || sym(1,3) != 13 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) != 12 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) != 15 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) != 13 || sym(3,2) != 14 || sym(3,3) != 11 || sym(3,4) != 19 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 15 || sym(4,3) != 19 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2 12 13 -1  8 )\n"
                                        "(  7 12 18 14 15 -2 )\n"
                                        "( -2 13 14 11 19  0 )\n"
                                        "(  5 -1 15 19  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 11 19 )
   {
      test_ = "Sparse matrix subtraction assignment test 3";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = -12;
         mat(0,1) = -11;
         mat(0,2) = -13;
         mat(0,3) = -15;
         mat(1,0) = -15;
         mat(1,1) = -14;
         mat(1,2) = -15;
         mat(1,3) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n( 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = -12;
         mat(0,1) = -15;
         mat(1,0) = -11;
         mat(1,1) = -14;
         mat(2,0) = -13;
         mat(2,1) = -15;
         mat(3,0) = -15;
         mat(3,1) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 11 ||
             sm(3,0) != 11 || sm(3,1) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 )\n"
                                        "( 18 14 )\n"
                                        "( 14 11 )\n"
                                        "( 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) != 12 || sym(2,5) != 13 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 18 || sym(3,5) != 14 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) != 12 || sym(4,3) != 18 || sym(4,4) != 14 || sym(4,5) != 11 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 13 || sym(5,3) != 14 || sym(5,4) != 11 || sym(5,5) != 19 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1 12 13 )\n"
                                        "( -2  0  1  5 18 14 )\n"
                                        "(  5 -1 12 18 14 11 )\n"
                                        "(  0  8 13 14 11 19 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 11 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix subtraction assignment test 4";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) =  -5;
         mat(0,1) = -18;
         mat(0,2) = -11;
         mat(0,3) = -10;
         mat(1,0) = -15;
         mat(1,1) = -14;
         mat(1,2) = -10;
         mat(1,3) = -14;
         mat(2,0) = -14;
         mat(2,1) = -12;
         mat(2,2) = -12;
         mat(2,3) =  -7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) != 11 ||
             sm(1,0) != 13 || sm(1,1) != 14 || sm(1,2) != 11 || sm(1,3) != 19 ||
             sm(2,0) != 19 || sm(2,1) != 11 || sm(2,2) != 12 || sm(2,3) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 18 14 11 )\n"
                                        "( 13 14 11 19 )\n"
                                        "( 19 11 12 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) =  -5;
         mat(0,1) = -15;
         mat(0,2) = -14;
         mat(1,0) = -18;
         mat(1,1) = -14;
         mat(1,2) = -12;
         mat(2,0) = -11;
         mat(2,1) = -10;
         mat(2,2) = -12;
         mat(3,0) = -10;
         mat(3,1) = -14;
         mat(3,2) =  -7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );
         sm -= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != 12 || sm(0,1) != 13 || sm(0,2) != 19 ||
             sm(1,0) != 18 || sm(1,1) != 14 || sm(1,2) != 11 ||
             sm(2,0) != 14 || sm(2,1) != 11 || sm(2,2) != 12 ||
             sm(3,0) != 11 || sm(3,1) != 19 || sm(3,2) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 12 13 19 )\n"
                                        "( 18 14 11 )\n"
                                        "( 14 11 12 )\n"
                                        "( 11 19 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 12 || sym(0,3) != 13 || sym(0,4) != 19 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) != 18 || sym(1,3) != 14 || sym(1,4) != 11 || sym(1,5) !=  8 ||
             sym(2,0) != 12 || sym(2,1) != 18 || sym(2,2) != 14 || sym(2,3) != 11 || sym(2,4) != 12 || sym(2,5) != -2 ||
             sym(3,0) != 13 || sym(3,1) != 14 || sym(3,2) != 11 || sym(3,3) != 19 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 19 || sym(4,1) != 11 || sym(4,2) != 12 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 12 13 19  0 )\n"
                                        "( -4  2 18 14 11  8 )\n"
                                        "( 12 18 14 11 12 -2 )\n"
                                        "( 13 14 11 19 14  0 )\n"
                                        "( 19 11 12 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      ( 12 18 14 15  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 22 17 11 19 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14 11  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 15 19  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix subtraction assignment test 5";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = -11;
         mat(0,1) = -22;
         mat(0,2) =  -7;
         mat(0,3) = -17;
         mat(1,0) = -26;
         mat(1,1) = -15;
         mat(1,2) = -11;
         mat(1,3) = -19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = -11;
         mat(0,1) = -22;
         mat(1,0) = -26;
         mat(1,1) = -15;
         mat(2,0) =  -7;
         mat(2,1) = -11;
         mat(3,0) = -17;
         mat(3,1) = -19;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 12 13 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7 12 18 14 15 -2 )
   // ( -2  0  1  5  7  0 )      ( -2 13 22 11 19  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 15 19  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix subtraction assignment test 6";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = -12;
         mat(0,1) = -15;
         mat(0,2) = -13;
         mat(0,3) = -15;
         mat(1,0) = -13;
         mat(1,1) = -21;
         mat(1,2) =  -6;
         mat(1,3) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = -12;
         mat(0,1) = -13;
         mat(1,0) = -15;
         mat(1,1) = -13;
         mat(2,0) = -21;
         mat(2,1) =  -6;
         mat(3,0) = -15;
         mat(3,1) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1 12 13 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 18 14 )
   // (  5 -1  0  7  1 -4 )      (  5 -1 12 18 14 11 )
   // (  0  8 -2  0 -4  7 )      (  0  8 13 14 22 19 )
   {
      test_ = "Sparse matrix subtraction assignment test 7";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = -12;
         mat(0,1) = -11;
         mat(0,2) = -13;
         mat(0,3) = -15;
         mat(1,0) = -15;
         mat(1,1) = -14;
         mat(1,2) = -26;
         mat(1,3) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = -12;
         mat(0,1) = -15;
         mat(1,0) = -11;
         mat(1,1) = -14;
         mat(2,0) = -13;
         mat(2,1) = -15;
         mat(3,0) = -26;
         mat(3,1) = -12;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 12 13 19  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2 18 14 11  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 12 18 14 11 12 -2 )
   // ( -2  0  1  5  7  0 )      ( 13 14 22 19 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 19 11 12 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix subtraction assignment test 8";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) =  -5;
         mat(0,1) = -18;
         mat(0,2) = -11;
         mat(0,3) = -10;
         mat(1,0) = -15;
         mat(1,1) = -14;
         mat(1,2) = -21;
         mat(1,3) = -14;
         mat(2,0) = -14;
         mat(2,1) = -12;
         mat(2,2) = -12;
         mat(2,3) =  -7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) =  -5;
         mat(0,1) = -15;
         mat(0,2) = -14;
         mat(1,0) = -18;
         mat(1,1) = -14;
         mat(1,2) = -12;
         mat(2,0) = -11;
         mat(2,1) = -10;
         mat(2,2) = -12;
         mat(3,0) = -21;
         mat(3,1) = -14;
         mat(3,2) =  -7;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Schur product assignment to a submatrix of a SymmetricMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment to a submatrix of a SymmetricMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename ST >  // Type of the symmetric matrix
void SubmatrixScalarTest::testSchurAssign()
{
   //=====================================================================================
   // Dense matrix Schur product assignment
   //=====================================================================================

   // (  1 -4  7 -2  5  0 )      ( 11 20 28 16  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 20 12  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 28  0  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 16  0  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix Schur product assignment test 1";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 11;
         mat(0,1) = -5;
         mat(0,2) =  4;
         mat(0,3) = -8;
         mat(1,0) = -5;
         mat(1,1) =  6;
         mat(1,2) = 99;
         mat(1,3) = 99;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) != 11 || sm(0,1) != 20 || sm(0,2) != 28 || sm(0,3) != 16 ||
             sm(1,0) != 20 || sm(1,1) != 12 || sm(1,2) !=  0 || sm(1,3) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 11 20 28 16 )\n( 20 12  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 11 || sym(0,1) != 20 || sym(0,2) != 28 || sym(0,3) != 16 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 20 || sym(1,1) != 12 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 28 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 16 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 11 20 28 16  5  0 )\n"
                                        "( 20 12  0  0 -1  8 )\n"
                                        "( 28  0  3  1  0 -2 )\n"
                                        "( 16  0  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 11;
         mat(0,1) = -5;
         mat(1,0) = -5;
         mat(1,1) =  6;
         mat(2,0) =  4;
         mat(2,1) = 99;
         mat(3,0) = -8;
         mat(3,1) = 99;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) != 11 || sm(0,1) != 20 ||
             sm(1,0) != 20 || sm(1,1) != 12 ||
             sm(2,0) != 28 || sm(2,1) !=  0 ||
             sm(3,0) != 16 || sm(3,1) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 11 20 )\n"
                                        "( 20 12 )\n"
                                        "( 28  0 )\n"
                                        "( 16  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 11 || sym(0,1) != 20 || sym(0,2) != 28 || sym(0,3) != 16 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 20 || sym(1,1) != 12 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 28 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 16 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 11 20 28 16  5  0 )\n"
                                        "( 20 12  0  0 -1  8 )\n"
                                        "( 28  0  3  1  0 -2 )\n"
                                        "( 16  0  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0 18 14  0 -2 )
   // ( -2  0  1  5  7  0 )      ( -2  0 14 20 21  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0 21  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix Schur product assignment test 2";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 99;
         mat(0,1) =  6;
         mat(0,2) = 14;
         mat(0,3) = 99;
         mat(1,0) = 99;
         mat(1,1) = 14;
         mat(1,2) =  4;
         mat(1,3) =  3;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) !=  0 ||
             sm(1,0) != 0 || sm(1,1) != 14 || sm(1,2) != 20 || sm(1,3) != 21 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 18 14  0 )\n( 0 14 20 21 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) != 14 || sym(3,3) != 20 || sym(3,4) != 21 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) != 21 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0 18 14  0 -2 )\n"
                                        "( -2  0 14 20 21  0 )\n"
                                        "(  5 -1  0 21  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 99;
         mat(0,1) = 99;
         mat(1,0) =  6;
         mat(1,1) = 14;
         mat(2,0) = 14;
         mat(2,1) =  4;
         mat(3,0) = 99;
         mat(3,1) =  3;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 20 ||
             sm(3,0) !=  0 || sm(3,1) != 21 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0  0 )\n"
                                        "( 18 14 )\n"
                                        "( 14 20 )\n"
                                        "(  0 21 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) != 14 || sym(3,3) != 20 || sym(3,4) != 21 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) != 21 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0 18 14  0 -2 )\n"
                                        "( -2  0 14 20 21  0 )\n"
                                        "(  5 -1  0 21  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1  0 16 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 21  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0 21 14 20 )
   // (  0  8 -2  0 -4  7 )      (  0  8 16  0 20 28 )
   {
      test_ = "Dense matrix Schur product assignment test 3";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 99;
         mat(0,1) =  3;
         mat(0,2) = 14;
         mat(0,3) = -5;
         mat(1,0) = -8;
         mat(1,1) = 99;
         mat(1,2) = -5;
         mat(1,3) =  4;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) !=  0 || sm(0,1) != 21 || sm(0,2) != 14 || sm(0,3) != 20 ||
             sm(1,0) != 16 || sm(1,1) !=  0 || sm(1,2) != 20 || sm(1,3) != 28 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0 21 14 20 )\n( 16  0 20 28 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != 16 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 21 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) != 21 || sym(4,4) != 14 || sym(4,5) != 20 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 16 || sym(5,3) !=  0 || sym(5,4) != 20 || sym(5,5) != 28 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1  0 16 )\n"
                                        "( -2  0  1  5 21  0 )\n"
                                        "(  5 -1  0 21 14 20 )\n"
                                        "(  0  8 16  0 20 28 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 99;
         mat(0,1) = -8;
         mat(1,0) =  3;
         mat(1,1) = 99;
         mat(2,0) = 14;
         mat(2,1) = -5;
         mat(3,0) = -5;
         mat(3,1) =  4;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) !=  0 || sm(0,1) != 16 ||
             sm(1,0) != 21 || sm(1,1) !=  0 ||
             sm(2,0) != 14 || sm(2,1) != 20 ||
             sm(3,0) != 20 || sm(3,1) != 28 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0 16 )\n"
                                        "( 21  0 )\n"
                                        "( 14 20 )\n"
                                        "( 20 28 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != 16 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 21 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) != 21 || sym(4,4) != 14 || sym(4,5) != 20 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 16 || sym(5,3) !=  0 || sym(5,4) != 20 || sym(5,5) != 28 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1  0 16 )\n"
                                        "( -2  0  1  5 21  0 )\n"
                                        "(  5 -1  0 21 14 20 )\n"
                                        "(  0  8 16  0 20 28 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 14 18 25  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0  7  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14  0 18 11  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 18  0 11 20 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 25  7  0 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix Schur product assignment test 4";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) =  2;
         mat(0,1) = 99;
         mat(0,2) =  6;
         mat(0,3) = 11;
         mat(1,0) = -9;
         mat(1,1) = 99;
         mat(1,2) = 11;
         mat(1,3) =  4;
         mat(2,0) =  5;
         mat(2,1) = -7;
         mat(2,2) = 99;
         mat(2,3) =  2;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) != 14 || sm(0,1) != 0 || sm(0,2) != 18 || sm(0,3) != 11 ||
             sm(1,0) != 18 || sm(1,1) != 0 || sm(1,2) != 11 || sm(1,3) != 20 ||
             sm(2,0) != 25 || sm(2,1) != 7 || sm(2,2) !=  0 || sm(2,3) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 14  0 18 11 )\n"
                                        "( 18  0 11 20 )\n"
                                        "( 25  7  0 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 14 || sym(0,3) != 18 || sym(0,4) != 25 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) !=  7 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) != 18 || sym(2,3) != 11 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 18 || sym(3,1) !=  0 || sym(3,2) != 11 || sym(3,3) != 20 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 25 || sym(4,1) !=  7 || sym(4,2) !=  0 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 14 18 25  0 )\n"
                                        "( -4  2  0  0  7  8 )\n"
                                        "( 14  0 18 11  0 -2 )\n"
                                        "( 18  0 11 20 14  0 )\n"
                                        "( 25  7  0 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) =  2;
         mat(0,1) = -9;
         mat(0,2) =  5;
         mat(1,0) = 99;
         mat(1,1) = 99;
         mat(1,2) = -7;
         mat(2,0) =  6;
         mat(2,1) = 11;
         mat(2,2) = 99;
         mat(3,0) = 11;
         mat(3,1) =  4;
         mat(3,2) =  2;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) != 14 || sm(0,1) != 18 || sm(0,2) != 25 ||
             sm(1,0) !=  0 || sm(1,1) !=  0 || sm(1,2) !=  7 ||
             sm(2,0) != 18 || sm(2,1) != 11 || sm(2,2) !=  0 ||
             sm(3,0) != 11 || sm(3,1) != 20 || sm(3,2) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 14 18 25 )\n"
                                        "(  0  0  7 )\n"
                                        "( 18 11  0 )\n"
                                        "( 11 20 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 14 || sym(0,3) != 18 || sym(0,4) != 25 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) !=  7 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) != 18 || sym(2,3) != 11 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 18 || sym(3,1) !=  0 || sym(3,2) != 11 || sym(3,3) != 20 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 25 || sym(4,1) !=  7 || sym(4,2) !=  0 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 14 18 25  0 )\n"
                                        "( -4  2  0  0  7  8 )\n"
                                        "( 14  0 18 11  0 -2 )\n"
                                        "( 18  0 11 20 14  0 )\n"
                                        "( 25  7  0 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      ( 11 20 28 16  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 24 12  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 28  0  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 16  0  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix Schur product assignment test 5";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 11;
         mat(0,1) = -5;
         mat(0,2) =  4;
         mat(0,3) = -8;
         mat(1,0) = -6;
         mat(1,1) =  6;
         mat(1,2) = 99;
         mat(1,3) = 99;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 11;
         mat(0,1) = -6;
         mat(1,0) = -5;
         mat(1,1) =  6;
         mat(2,0) =  4;
         mat(2,1) = 99;
         mat(3,0) = -8;
         mat(3,1) = 99;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0 18 14  0 -2 )
   // ( -2  0  1  5  7  0 )      ( -2  0 22 20 21  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0 21  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix Schur product assignment test 6";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 99;
         mat(0,1) =  6;
         mat(0,2) = 14;
         mat(0,3) = 99;
         mat(1,0) = 99;
         mat(1,1) = 22;
         mat(1,2) =  4;
         mat(1,3) =  3;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 99;
         mat(0,1) = 99;
         mat(1,0) =  6;
         mat(1,1) = 22;
         mat(2,0) = 14;
         mat(2,1) =  4;
         mat(3,0) = 99;
         mat(3,1) =  3;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1  0 16 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 21  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0 21 14 20 )
   // (  0  8 -2  0 -4  7 )      (  0  8 16  0 24 28 )
   {
      test_ = "Dense matrix Schur product assignment test 7";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = 99;
         mat(0,1) =  3;
         mat(0,2) = 14;
         mat(0,3) = -5;
         mat(1,0) = -8;
         mat(1,1) = 99;
         mat(1,2) = -6;
         mat(1,3) =  4;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = 99;
         mat(0,1) = -8;
         mat(1,0) =  3;
         mat(1,1) = 99;
         mat(2,0) = 14;
         mat(2,1) = -6;
         mat(3,0) = -5;
         mat(3,1) =  4;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 14 18 25  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0  7  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14  0 18 11  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 18  0 22 20 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 25  7  0 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Dense matrix Schur product assignment test 8";

      {
         blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) =  2;
         mat(0,1) = 99;
         mat(0,2) =  6;
         mat(0,3) = 11;
         mat(1,0) = -9;
         mat(1,1) = 99;
         mat(1,2) = 22;
         mat(1,3) =  4;
         mat(2,0) =  5;
         mat(2,1) = -7;
         mat(2,2) = 99;
         mat(2,3) =  2;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<int,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) =  2;
         mat(0,1) = -9;
         mat(0,2) =  5;
         mat(1,0) = 99;
         mat(1,1) = 99;
         mat(1,2) = -7;
         mat(2,0) =  6;
         mat(2,1) = 22;
         mat(2,2) = 99;
         mat(3,0) = 11;
         mat(3,1) =  4;
         mat(3,2) =  2;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Sparse matrix Schur product assignment
   //=====================================================================================

   // (  1 -4  7 -2  5  0 )      ( 11 20 28 16  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 20 12  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 28  0  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 16  0  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix Schur product assignment test 1";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 11;
         mat(0,1) = -5;
         mat(0,2) =  4;
         mat(0,3) = -8;
         mat(1,0) = -5;
         mat(1,1) =  6;
         mat(1,2) = 99;
         mat(1,3) = 99;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) != 11 || sm(0,1) != 20 || sm(0,2) != 28 || sm(0,3) != 16 ||
             sm(1,0) != 20 || sm(1,1) != 12 || sm(1,2) !=  0 || sm(1,3) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 11 20 28 16 )\n( 20 12  0  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 11 || sym(0,1) != 20 || sym(0,2) != 28 || sym(0,3) != 16 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 20 || sym(1,1) != 12 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 28 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 16 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 11 20 28 16  5  0 )\n"
                                        "( 20 12  0  0 -1  8 )\n"
                                        "( 28  0  3  1  0 -2 )\n"
                                        "( 16  0  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 11;
         mat(0,1) = -5;
         mat(1,0) = -5;
         mat(1,1) =  6;
         mat(2,0) =  4;
         mat(2,1) = 99;
         mat(3,0) = -8;
         mat(3,1) = 99;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) != 11 || sm(0,1) != 20 ||
             sm(1,0) != 20 || sm(1,1) != 12 ||
             sm(2,0) != 28 || sm(2,1) !=  0 ||
             sm(3,0) != 16 || sm(3,1) !=  0 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 11 20 )\n"
                                        "( 20 12 )\n"
                                        "( 28  0 )\n"
                                        "( 16  0 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != 11 || sym(0,1) != 20 || sym(0,2) != 28 || sym(0,3) != 16 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != 20 || sym(1,1) != 12 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) != 28 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 16 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) !=  7 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) !=  7 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( 11 20 28 16  5  0 )\n"
                                        "( 20 12  0  0 -1  8 )\n"
                                        "( 28  0  3  1  0 -2 )\n"
                                        "( 16  0  1  5  7  0 )\n"
                                        "(  5 -1  0  7  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0 18 14  0 -2 )
   // ( -2  0  1  5  7  0 )      ( -2  0 14 20 21  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0 21  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix Schur product assignment test 2";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 99;
         mat(0,1) =  6;
         mat(0,2) = 14;
         mat(0,3) = 99;
         mat(1,0) = 99;
         mat(1,1) = 14;
         mat(1,2) =  4;
         mat(1,3) =  3;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) != 0 || sm(0,1) != 18 || sm(0,2) != 14 || sm(0,3) !=  0 ||
             sm(1,0) != 0 || sm(1,1) != 14 || sm(1,2) != 20 || sm(1,3) != 21 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 0 18 14  0 )\n( 0 14 20 21 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) != 14 || sym(3,3) != 20 || sym(3,4) != 21 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) != 21 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0 18 14  0 -2 )\n"
                                        "( -2  0 14 20 21  0 )\n"
                                        "(  5 -1  0 21  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 99;
         mat(0,1) = 99;
         mat(1,0) =  6;
         mat(1,1) = 14;
         mat(2,0) = 14;
         mat(2,1) =  4;
         mat(3,0) = 99;
         mat(3,1) =  3;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) !=  0 || sm(0,1) !=  0 ||
             sm(1,0) != 18 || sm(1,1) != 14 ||
             sm(2,0) != 14 || sm(2,1) != 20 ||
             sm(3,0) !=  0 || sm(3,1) != 21 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0  0 )\n"
                                        "( 18 14 )\n"
                                        "( 14 20 )\n"
                                        "(  0 21 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) != 18 || sym(2,3) != 14 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) != 14 || sym(3,3) != 20 || sym(3,4) != 21 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) != 21 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0 18 14  0 -2 )\n"
                                        "( -2  0 14 20 21  0 )\n"
                                        "(  5 -1  0 21  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1  0 16 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 21  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0 21 14 20 )
   // (  0  8 -2  0 -4  7 )      (  0  8 16  0 20 28 )
   {
      test_ = "Sparse matrix Schur product assignment test 3";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 99;
         mat(0,1) =  3;
         mat(0,2) = 14;
         mat(0,3) = -5;
         mat(1,0) = -8;
         mat(1,1) = 99;
         mat(1,2) = -5;
         mat(1,3) =  4;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) !=  0 || sm(0,1) != 21 || sm(0,2) != 14 || sm(0,3) != 20 ||
             sm(1,0) != 16 || sm(1,1) !=  0 || sm(1,2) != 20 || sm(1,3) != 28 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0 21 14 20 )\n( 16  0 20 28 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != 16 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 21 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) != 21 || sym(4,4) != 14 || sym(4,5) != 20 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 16 || sym(5,3) !=  0 || sym(5,4) != 20 || sym(5,5) != 28 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1  0 16 )\n"
                                        "( -2  0  1  5 21  0 )\n"
                                        "(  5 -1  0 21 14 20 )\n"
                                        "(  0  8 16  0 20 28 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 99;
         mat(0,1) = -8;
         mat(1,0) =  3;
         mat(1,1) = 99;
         mat(2,0) = 14;
         mat(2,1) = -5;
         mat(3,0) = -5;
         mat(3,1) =  4;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) !=  0 || sm(0,1) != 16 ||
             sm(1,0) != 21 || sm(1,1) !=  0 ||
             sm(2,0) != 14 || sm(2,1) != 20 ||
             sm(3,0) != 20 || sm(3,1) != 28 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n(  0 16 )\n"
                                        "( 21  0 )\n"
                                        "( 14 20 )\n"
                                        "( 20 28 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) !=  7 || sym(0,3) != -2 || sym(0,4) !=  5 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) != -1 || sym(1,5) !=  8 ||
             sym(2,0) !=  7 || sym(2,1) !=  0 || sym(2,2) !=  3 || sym(2,3) !=  1 || sym(2,4) !=  0 || sym(2,5) != 16 ||
             sym(3,0) != -2 || sym(3,1) !=  0 || sym(3,2) !=  1 || sym(3,3) !=  5 || sym(3,4) != 21 || sym(3,5) !=  0 ||
             sym(4,0) !=  5 || sym(4,1) != -1 || sym(4,2) !=  0 || sym(4,3) != 21 || sym(4,4) != 14 || sym(4,5) != 20 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != 16 || sym(5,3) !=  0 || sym(5,4) != 20 || sym(5,5) != 28 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4  7 -2  5  0 )\n"
                                        "( -4  2  0  0 -1  8 )\n"
                                        "(  7  0  3  1  0 16 )\n"
                                        "( -2  0  1  5 21  0 )\n"
                                        "(  5 -1  0 21 14 20 )\n"
                                        "(  0  8 16  0 20 28 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 14 18 25  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0  7  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14  0 18 11  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 18  0 11 20 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 25  7  0 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix Schur product assignment test 4";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) =  2;
         mat(0,1) = 99;
         mat(0,2) =  6;
         mat(0,3) = 11;
         mat(1,0) = -9;
         mat(1,1) = 99;
         mat(1,2) = 11;
         mat(1,3) =  4;
         mat(2,0) =  5;
         mat(2,1) = -7;
         mat(2,2) = 99;
         mat(2,3) =  2;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) != 14 || sm(0,1) != 0 || sm(0,2) != 18 || sm(0,3) != 11 ||
             sm(1,0) != 18 || sm(1,1) != 0 || sm(1,2) != 11 || sm(1,3) != 20 ||
             sm(2,0) != 25 || sm(2,1) != 7 || sm(2,2) !=  0 || sm(2,3) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 14  0 18 11 )\n"
                                        "( 18  0 11 20 )\n"
                                        "( 25  7  0 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 14 || sym(0,3) != 18 || sym(0,4) != 25 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) !=  7 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) != 18 || sym(2,3) != 11 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 18 || sym(3,1) !=  0 || sym(3,2) != 11 || sym(3,3) != 20 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 25 || sym(4,1) !=  7 || sym(4,2) !=  0 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 14 18 25  0 )\n"
                                        "( -4  2  0  0  7  8 )\n"
                                        "( 14  0 18 11  0 -2 )\n"
                                        "( 18  0 11 20 14  0 )\n"
                                        "( 25  7  0 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) =  2;
         mat(0,1) = -9;
         mat(0,2) =  5;
         mat(1,0) = 99;
         mat(1,1) = 99;
         mat(1,2) = -7;
         mat(2,0) =  6;
         mat(2,1) = 11;
         mat(2,2) = 99;
         mat(3,0) = 11;
         mat(3,1) =  4;
         mat(3,2) =  2;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );
         sm %= mat;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 26UL );

         if( sm(0,0) != 14 || sm(0,1) != 18 || sm(0,2) != 25 ||
             sm(1,0) !=  0 || sm(1,1) !=  0 || sm(1,2) !=  7 ||
             sm(2,0) != 18 || sm(2,1) != 11 || sm(2,2) !=  0 ||
             sm(3,0) != 11 || sm(3,1) != 20 || sm(3,2) != 14 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( 14 18 25 )\n"
                                        "(  0  0  7 )\n"
                                        "( 18 11  0 )\n"
                                        "( 11 20 14 )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) !=  1 || sym(0,1) != -4 || sym(0,2) != 14 || sym(0,3) != 18 || sym(0,4) != 25 || sym(0,5) !=  0 ||
             sym(1,0) != -4 || sym(1,1) !=  2 || sym(1,2) !=  0 || sym(1,3) !=  0 || sym(1,4) !=  7 || sym(1,5) !=  8 ||
             sym(2,0) != 14 || sym(2,1) !=  0 || sym(2,2) != 18 || sym(2,3) != 11 || sym(2,4) !=  0 || sym(2,5) != -2 ||
             sym(3,0) != 18 || sym(3,1) !=  0 || sym(3,2) != 11 || sym(3,3) != 20 || sym(3,4) != 14 || sym(3,5) !=  0 ||
             sym(4,0) != 25 || sym(4,1) !=  7 || sym(4,2) !=  0 || sym(4,3) != 14 || sym(4,4) !=  1 || sym(4,5) != -4 ||
             sym(5,0) !=  0 || sym(5,1) !=  8 || sym(5,2) != -2 || sym(5,3) !=  0 || sym(5,4) != -4 || sym(5,5) !=  7 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n(  1 -4 14 18 25  0 )\n"
                                        "( -4  2  0  0  7  8 )\n"
                                        "( 14  0 18 11  0 -2 )\n"
                                        "( 18  0 11 20 14  0 )\n"
                                        "( 25  7  0 14  1 -4 )\n"
                                        "(  0  8 -2  0 -4  7 )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // (  1 -4  7 -2  5  0 )      ( 11 20 28 16  5  0 )
   // ( -4  2  0  0 -1  8 )      ( 24 12  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 28  0  3  1  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 16  0  1  5  7  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0  7  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix Schur product assignment test 5";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 11;
         mat(0,1) = -5;
         mat(0,2) =  4;
         mat(0,3) = -8;
         mat(1,0) = -6;
         mat(1,1) =  6;
         mat(1,2) = 99;
         mat(1,3) = 99;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 11;
         mat(0,1) = -6;
         mat(1,0) = -5;
         mat(1,1) =  6;
         mat(2,0) =  4;
         mat(2,1) = 99;
         mat(3,0) = -8;
         mat(3,1) = 99;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0 18 14  0 -2 )
   // ( -2  0  1  5  7  0 )      ( -2  0 22 20 21  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0 21  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix Schur product assignment test 6";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 99;
         mat(0,1) =  6;
         mat(0,2) = 14;
         mat(0,3) = 99;
         mat(1,0) = 99;
         mat(1,1) = 22;
         mat(1,2) =  4;
         mat(1,3) =  3;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 99;
         mat(0,1) = 99;
         mat(1,0) =  6;
         mat(1,1) = 22;
         mat(2,0) = 14;
         mat(2,1) =  4;
         mat(3,0) = 99;
         mat(3,1) =  3;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4  7 -2  5  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0 -1  8 )
   // (  7  0  3  1  0 -2 )  =>  (  7  0  3  1  0 16 )
   // ( -2  0  1  5  7  0 )      ( -2  0  1  5 21  0 )
   // (  5 -1  0  7  1 -4 )      (  5 -1  0 21 14 20 )
   // (  0  8 -2  0 -4  7 )      (  0  8 16  0 24 28 )
   {
      test_ = "Sparse matrix Schur product assignment test 7";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = 99;
         mat(0,1) =  3;
         mat(0,2) = 14;
         mat(0,3) = -5;
         mat(1,0) = -8;
         mat(1,1) = 99;
         mat(1,2) = -6;
         mat(1,3) =  4;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = 99;
         mat(0,1) = -8;
         mat(1,0) =  3;
         mat(1,1) = 99;
         mat(2,0) = 14;
         mat(2,1) = -6;
         mat(3,0) = -5;
         mat(3,1) =  4;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // (  1 -4  7 -2  5  0 )      (  1 -4 14 18 25  0 )
   // ( -4  2  0  0 -1  8 )      ( -4  2  0  0  7  8 )
   // (  7  0  3  1  0 -2 )  =>  ( 14  0 18 11  0 -2 )
   // ( -2  0  1  5  7  0 )      ( 18  0 22 20 14  0 )
   // (  5 -1  0  7  1 -4 )      ( 25  7  0 14  1 -4 )
   // (  0  8 -2  0 -4  7 )      (  0  8 -2  0 -4  7 )
   {
      test_ = "Sparse matrix Schur product assignment test 8";

      {
         blaze::CompressedMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) =  2;
         mat(0,1) = 99;
         mat(0,2) =  6;
         mat(0,3) = 11;
         mat(1,0) = -9;
         mat(1,1) = 99;
         mat(1,2) = 22;
         mat(1,3) =  4;
         mat(2,0) =  5;
         mat(2,1) = -7;
         mat(2,2) = 99;
         mat(2,3) =  2;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<int,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) =  2;
         mat(0,1) = -9;
         mat(0,2) =  5;
         mat(1,0) = 99;
         mat(1,1) = 99;
         mat(1,2) = -7;
         mat(2,0) =  6;
         mat(2,1) = 22;
         mat(2,2) = 99;
         mat(3,0) = 11;
         mat(3,1) =  4;
         mat(3,2) =  2;

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
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
void SubmatrixScalarTest::checkRows( const Type& matrix, size_t expectedRows ) const
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
void SubmatrixScalarTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
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
void SubmatrixScalarTest::checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const
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
/*!\brief Initializing the given symmetric matrix.
//
// \return void
//
// This function is called before each test case to initialize the given symmetric matrix.
*/
template< typename ST >
void SubmatrixScalarTest::init( ST& sym )
{
   sym.resize( 6UL );
   sym(0,0) =  1;
   sym(0,1) = -4;
   sym(0,2) =  7;
   sym(0,3) = -2;
   sym(0,4) =  5;
   sym(1,1) =  2;
   sym(1,4) = -1;
   sym(1,5) =  8;
   sym(2,2) =  3;
   sym(2,3) =  1;
   sym(2,5) = -2;
   sym(3,3) =  5;
   sym(3,4) =  7;
   sym(4,4) =  1;
   sym(4,5) = -4;
   sym(5,5) =  7;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the assignment to a submatrix of a scalar SymmetricMatrix.
//
// \return void
*/
void runTest()
{
   SubmatrixScalarTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the SymmetricMatrix submatrix scalar test.
*/
#define RUN_SYMMETRICMATRIX_SUBMATRIXSCALAR_TEST \
   blazetest::mathtest::adaptors::symmetricmatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace symmetricmatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
