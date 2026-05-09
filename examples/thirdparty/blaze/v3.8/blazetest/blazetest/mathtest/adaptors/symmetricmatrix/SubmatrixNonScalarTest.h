//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/symmetricmatrix/SubmatrixNonScalarTest.h
//  \brief Header file for the SymmetricMatrix submatrix non-scalar test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_SYMMETRICMATRIX_SUBMATRIXNONSCALARTEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_SYMMETRICMATRIX_SUBMATRIXNONSCALARTEST_H_


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
/*!\brief Auxiliary class for assignment tests to a submatrix of a non-scalar SymmetricMatrix.
//
// This class performs assignment tests to a submatrix of a SymmetricMatrix with non-scalar
// element type. It performs a series of both compile time as well as runtime tests.
*/
class SubmatrixNonScalarTest
{
 private:
   //**Type definitions****************************************************************************
   //! Type of a resizable, non-scalar element.
   using VT = blaze::DynamicVector<int,blaze::rowVector>;

   //! Type of the dense non-scalar symmetric matrix.
   using DST = blaze::SymmetricMatrix< blaze::DynamicMatrix<VT,blaze::rowMajor> >;

   //! Opposite dense non-scalar symmetric matrix type.
   using DOST = DST::OppositeType;

   //! Type of the sparse non-scalar symmetric matrix.
   using SST = blaze::SymmetricMatrix< blaze::CompressedMatrix<VT,blaze::rowMajor> >;

   //! Opposite sparse non-scalar symmetric matrix type.
   using SOST = SST::OppositeType;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit SubmatrixNonScalarTest();
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
   template< typename ST > void testAssignment();

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
   inline VT vec( int value );
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
void SubmatrixNonScalarTest::testAssignment()
{
   //=====================================================================================
   // Dense matrix assignment
   //=====================================================================================

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( ( 12 ) ( 18 ) ( 14 ) ( 15 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( 18 ) ( 17 ) ( 11 ) ( 19 ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( ( 14 ) ( 11 ) (  3 ) (  1 ) (    ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( 15 ) ( 19 ) (  1 ) (  5 ) (  7 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Dense matrix assignment test 1";

      {
         blaze::DynamicMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 15 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 17 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 18 ) || sm(0,2) != vec( 14 ) || sm(0,3) != vec( 15 ) ||
             sm(1,0) != vec( 18 ) || sm(1,1) != vec( 17 ) || sm(1,2) != vec( 11 ) || sm(1,3) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 15 ) )\n"
                                        "( ( 18 ) ( 17 ) ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec( 12 )  || sym(0,1) != vec( 18 ) || sym(0,2) != vec( 14 )  || sym(0,3) != vec( 15 )  || sym(0,4) != vec(  5 )  || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( 18 )  || sym(1,1) != vec( 17 ) || sym(1,2) != vec( 11 )  || sym(1,3) != vec( 19 )  || sym(1,4) != vec( -1 )  || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec( 14 )  || sym(2,1) != vec( 11 ) || sym(2,2) != vec(  3 )  || sym(2,3) != vec(  1 )  || !isDefault( sym(2,4) ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( 15 )  || sym(3,1) != vec( 19 ) || sym(3,2) != vec(  1 )  || sym(3,3) != vec(  5 )  || sym(3,4) != vec(  7 )  || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 ) || !isDefault( sym(4,2) ) || sym(4,3) != vec(  7 )  || sym(4,4) != vec(  1 )  || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 )  || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 )  || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 15 ) (  5 ) (    ) )\n"
                                        "( ( 18 ) ( 17 ) ( 11 ) ( 19 ) ( -1 ) (  8 ) )\n"
                                        "( ( 14 ) ( 11 ) (  3 ) (  1 ) (    ) ( -2 ) )\n"
                                        "( ( 15 ) ( 19 ) (  1 ) (  5 ) (  7 ) (    ) )\n"
                                        "( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 17 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 15 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 18 ) ||
             sm(1,0) != vec( 18 ) || sm(1,1) != vec( 17 ) ||
             sm(2,0) != vec( 14 ) || sm(2,1) != vec( 11 ) ||
             sm(3,0) != vec( 15 ) || sm(3,1) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) )\n"
                                        "( ( 18 ) ( 17 ) )\n"
                                        "( ( 14 ) ( 11 ) )\n"
                                        "( ( 15 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec( 12 )  || sym(0,1) != vec( 18 ) || sym(0,2) != vec( 14 )  || sym(0,3) != vec( 15 )  || sym(0,4) != vec(  5 )  || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( 18 )  || sym(1,1) != vec( 17 ) || sym(1,2) != vec( 11 )  || sym(1,3) != vec( 19 )  || sym(1,4) != vec( -1 )  || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec( 14 )  || sym(2,1) != vec( 11 ) || sym(2,2) != vec(  3 )  || sym(2,3) != vec(  1 )  || !isDefault( sym(2,4) ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( 15 )  || sym(3,1) != vec( 19 ) || sym(3,2) != vec(  1 )  || sym(3,3) != vec(  5 )  || sym(3,4) != vec(  7 )  || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 ) || !isDefault( sym(4,2) ) || sym(4,3) != vec(  7 )  || sym(4,4) != vec(  1 )  || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 )  || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 )  || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 15 ) (  5 ) (    ) )\n"
                                        "( ( 18 ) ( 17 ) ( 11 ) ( 19 ) ( -1 ) (  8 ) )\n"
                                        "( ( 14 ) ( 11 ) (  3 ) (  1 ) (    ) ( -2 ) )\n"
                                        "( ( 15 ) ( 19 ) (  1 ) (  5 ) (  7 ) (    ) )\n"
                                        "( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) ( 12 ) ( 13 ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( (  7 ) ( 12 ) ( 18 ) ( 14 ) ( 15 ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( -2 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) ( 15 ) ( 19 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Dense matrix assignment test 2";

      {
         blaze::DynamicMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 15 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 18 ) || sm(0,2) != vec( 14 ) || sm(0,3) != vec( 15 ) ||
             sm(1,0) != vec( 13 ) || sm(1,1) != vec( 14 ) || sm(1,2) != vec( 11 ) || sm(1,3) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 15 ) )\n"
                                        "( ( 13 ) ( 14 ) ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 ) || sym(0,2) != vec(  7 ) || sym(0,3) != vec( -2 )  || sym(0,4) != vec(  5 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 ) || sym(1,2) != vec( 12 ) || sym(1,3) != vec( 13 )  || sym(1,4) != vec( -1 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec(  7 )  || sym(2,1) != vec( 12 ) || sym(2,2) != vec( 18 ) || sym(2,3) != vec( 14 )  || sym(2,4) != vec( 15 ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( -2 )  || sym(3,1) != vec( 13 ) || sym(3,2) != vec( 14 ) || sym(3,3) != vec( 11 )  || sym(3,4) != vec( 19 ) || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 ) || sym(4,2) != vec( 15 ) || sym(4,3) != vec( 19 )  || sym(4,4) != vec(  1 ) || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 ) || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 ) || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) ( 12 ) ( 13 ) ( -1 ) (  8 ) )\n"
                                        "( (  7 ) ( 12 ) ( 18 ) ( 14 ) ( 15 ) ( -2 ) )\n"
                                        "( ( -2 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) (    ) )\n"
                                        "( (  5 ) ( -1 ) ( 15 ) ( 19 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 15 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 13 ) ||
             sm(1,0) != vec( 18 ) || sm(1,1) != vec( 14 ) ||
             sm(2,0) != vec( 14 ) || sm(2,1) != vec( 11 ) ||
             sm(3,0) != vec( 15 ) || sm(3,1) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 13 ) )\n"
                                        "( ( 18 ) ( 14 ) )\n"
                                        "( ( 14 ) ( 11 ) )\n"
                                        "( ( 15 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 ) || sym(0,2) != vec(  7 ) || sym(0,3) != vec( -2 )  || sym(0,4) != vec(  5 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 ) || sym(1,2) != vec( 12 ) || sym(1,3) != vec( 13 )  || sym(1,4) != vec( -1 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec(  7 )  || sym(2,1) != vec( 12 ) || sym(2,2) != vec( 18 ) || sym(2,3) != vec( 14 )  || sym(2,4) != vec( 15 ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( -2 )  || sym(3,1) != vec( 13 ) || sym(3,2) != vec( 14 ) || sym(3,3) != vec( 11 )  || sym(3,4) != vec( 19 ) || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 ) || sym(4,2) != vec( 15 ) || sym(4,3) != vec( 19 )  || sym(4,4) != vec(  1 ) || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 ) || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 ) || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) ( 12 ) ( 13 ) ( -1 ) (  8 ) )\n"
                                        "( (  7 ) ( 12 ) ( 18 ) ( 14 ) ( 15 ) ( -2 ) )\n"
                                        "( ( -2 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) (    ) )\n"
                                        "( (  5 ) ( -1 ) ( 15 ) ( 19 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( (  7 ) (    ) (  3 ) (  1 ) ( 12 ) ( 13 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( -2 ) (    ) (  1 ) (  5 ) ( 18 ) ( 14 ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) ( 12 ) ( 18 ) ( 14 ) ( 11 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) )
   {
      test_ = "Dense matrix assignment test 3";

      {
         blaze::DynamicMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 11 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 18 ) || sm(0,2) != vec( 14 ) || sm(0,3) != vec( 11 ) ||
             sm(1,0) != vec( 13 ) || sm(1,1) != vec( 14 ) || sm(1,2) != vec( 11 ) || sm(1,3) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 11 ) )\n"
                                        "( ( 13 ) ( 14 ) ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 )  || sym(0,2) != vec(  7 )  || sym(0,3) != vec( -2 )  || sym(0,4) != vec(  5 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) || !isDefault( sym(1,3) ) || sym(1,4) != vec( -1 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec(  7 )  || !isDefault( sym(2,1) ) || sym(2,2) != vec(  3 )  || sym(2,3) != vec(  1 )  || sym(2,4) != vec( 12 ) || sym(2,5) != vec( 13 )  ||
             sym(3,0) != vec( -2 )  || !isDefault( sym(3,1) ) || sym(3,2) != vec(  1 )  || sym(3,3) != vec(  5 )  || sym(3,4) != vec( 18 ) || sym(3,5) != vec( 14 )  ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 )  || sym(4,2) != vec( 12 )  || sym(4,3) != vec( 18 )  || sym(4,4) != vec( 14 ) || sym(4,5) != vec( 11 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 )  || sym(5,2) != vec( 13 )  || sym(5,3) != vec( 14 )  || sym(5,4) != vec( 11 ) || sym(5,5) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )\n"
                                        "( (  7 ) (    ) (  3 ) (  1 ) ( 12 ) ( 13 ) )\n"
                                        "( ( -2 ) (    ) (  1 ) (  5 ) ( 18 ) ( 14 ) )\n"
                                        "( (  5 ) ( -1 ) ( 12 ) ( 18 ) ( 14 ) ( 11 ) )\n"
                                        "( (    ) (  8 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 11 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 13 ) ||
             sm(1,0) != vec( 18 ) || sm(1,1) != vec( 14 ) ||
             sm(2,0) != vec( 14 ) || sm(2,1) != vec( 11 ) ||
             sm(3,0) != vec( 11 ) || sm(3,1) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 13 ) )\n"
                                        "( ( 18 ) ( 14 ) )\n"
                                        "( ( 14 ) ( 11 ) )\n"
                                        "( ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 )  || sym(0,2) != vec(  7 )  || sym(0,3) != vec( -2 )  || sym(0,4) != vec(  5 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) || !isDefault( sym(1,3) ) || sym(1,4) != vec( -1 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec(  7 )  || !isDefault( sym(2,1) ) || sym(2,2) != vec(  3 )  || sym(2,3) != vec(  1 )  || sym(2,4) != vec( 12 ) || sym(2,5) != vec( 13 )  ||
             sym(3,0) != vec( -2 )  || !isDefault( sym(3,1) ) || sym(3,2) != vec(  1 )  || sym(3,3) != vec(  5 )  || sym(3,4) != vec( 18 ) || sym(3,5) != vec( 14 )  ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 )  || sym(4,2) != vec( 12 )  || sym(4,3) != vec( 18 )  || sym(4,4) != vec( 14 ) || sym(4,5) != vec( 11 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 )  || sym(5,2) != vec( 13 )  || sym(5,3) != vec( 14 )  || sym(5,4) != vec( 11 ) || sym(5,5) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )\n"
                                        "( (  7 ) (    ) (  3 ) (  1 ) ( 12 ) ( 13 ) )\n"
                                        "( ( -2 ) (    ) (  1 ) (  5 ) ( 18 ) ( 14 ) )\n"
                                        "( (  5 ) ( -1 ) ( 12 ) ( 18 ) ( 14 ) ( 11 ) )\n"
                                        "( (    ) (  8 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) ( 12 ) ( 13 ) ( 19 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) ( 18 ) ( 14 ) ( 11 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( ( 12 ) ( 18 ) ( 14 ) ( 11 ) ( 12 ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( 13 ) ( 14 ) ( 11 ) ( 19 ) ( 14 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( ( 19 ) ( 11 ) ( 12 ) ( 14 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Dense matrix assignment test 4";

      {
         blaze::DynamicMatrix<VT,blaze::rowMajor> tmp( 3UL, 4UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 11 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );
         tmp(2,0) = vec( 19 );
         tmp(2,1) = vec( 11 );
         tmp(2,2) = vec( 12 );
         tmp(2,3) = vec( 14 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 18 ) || sm(0,2) != vec( 14 ) || sm(0,3) != vec( 11 ) ||
             sm(1,0) != vec( 13 ) || sm(1,1) != vec( 14 ) || sm(1,2) != vec( 11 ) || sm(1,3) != vec( 19 ) ||
             sm(2,0) != vec( 19 ) || sm(2,1) != vec( 11 ) || sm(2,2) != vec( 12 ) || sm(2,3) != vec( 14 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 11 ) )\n"
                                        "( ( 13 ) ( 14 ) ( 11 ) ( 19 ) )\n"
                                        "( ( 19 ) ( 11 ) ( 12 ) ( 14 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 ) || sym(0,2) != vec( 12 ) || sym(0,3) != vec( 13 )  || sym(0,4) != vec( 19 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 ) || sym(1,2) != vec( 18 ) || sym(1,3) != vec( 14 )  || sym(1,4) != vec( 11 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec( 12 )  || sym(2,1) != vec( 18 ) || sym(2,2) != vec( 14 ) || sym(2,3) != vec( 11 )  || sym(2,4) != vec( 12 ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( 13 )  || sym(3,1) != vec( 14 ) || sym(3,2) != vec( 11 ) || sym(3,3) != vec( 19 )  || sym(3,4) != vec( 14 ) || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec( 19 )  || sym(4,1) != vec( 11 ) || sym(4,2) != vec( 12 ) || sym(4,3) != vec( 14 )  || sym(4,4) != vec(  1 ) || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 ) || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 ) || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) ( 12 ) ( 13 ) ( 19 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) ( 18 ) ( 14 ) ( 11 ) (  8 ) )\n"
                                        "( ( 12 ) ( 18 ) ( 14 ) ( 11 ) ( 12 ) ( -2 ) )\n"
                                        "( ( 13 ) ( 14 ) ( 11 ) ( 19 ) ( 14 ) (    ) )\n"
                                        "( ( 19 ) ( 11 ) ( 12 ) ( 14 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<VT,blaze::rowMajor> tmp( 4UL, 3UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(0,2) = vec( 19 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 11 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(2,2) = vec( 12 );
         tmp(3,0) = vec( 11 );
         tmp(3,1) = vec( 19 );
         tmp(3,2) = vec( 14 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 13 ) || sm(0,2) != vec( 19 ) ||
             sm(1,0) != vec( 18 ) || sm(1,1) != vec( 14 ) || sm(1,2) != vec( 11 ) ||
             sm(2,0) != vec( 14 ) || sm(2,1) != vec( 11 ) || sm(2,2) != vec( 12 ) ||
             sm(3,0) != vec( 11 ) || sm(3,1) != vec( 19 ) || sm(3,2) != vec( 14 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 13 ) ( 19 ) )\n"
                                        "( ( 18 ) ( 14 ) ( 11 ) )\n"
                                        "( ( 14 ) ( 11 ) ( 12 ) )\n"
                                        "( ( 11 ) ( 19 ) ( 14 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 ) || sym(0,2) != vec( 12 ) || sym(0,3) != vec( 13 )  || sym(0,4) != vec( 19 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 ) || sym(1,2) != vec( 18 ) || sym(1,3) != vec( 14 )  || sym(1,4) != vec( 11 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec( 12 )  || sym(2,1) != vec( 18 ) || sym(2,2) != vec( 14 ) || sym(2,3) != vec( 11 )  || sym(2,4) != vec( 12 ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( 13 )  || sym(3,1) != vec( 14 ) || sym(3,2) != vec( 11 ) || sym(3,3) != vec( 19 )  || sym(3,4) != vec( 14 ) || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec( 19 )  || sym(4,1) != vec( 11 ) || sym(4,2) != vec( 12 ) || sym(4,3) != vec( 14 )  || sym(4,4) != vec(  1 ) || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 ) || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 ) || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) ( 12 ) ( 13 ) ( 19 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) ( 18 ) ( 14 ) ( 11 ) (  8 ) )\n"
                                        "( ( 12 ) ( 18 ) ( 14 ) ( 11 ) ( 12 ) ( -2 ) )\n"
                                        "( ( 13 ) ( 14 ) ( 11 ) ( 19 ) ( 14 ) (    ) )\n"
                                        "( ( 19 ) ( 11 ) ( 12 ) ( 14 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( ( 12 ) ( 18 ) ( 14 ) ( 15 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( 22 ) ( 17 ) ( 11 ) ( 19 ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( ( 14 ) ( 11 ) (  3 ) (  1 ) (    ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( 15 ) ( 19 ) (  1 ) (  5 ) (  7 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Dense matrix assignment test 5";

      {
         blaze::DynamicMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 15 );
         tmp(1,0) = vec( 22 );
         tmp(1,1) = vec( 17 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );

         try {
            sm = tmp;

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
         blaze::DynamicMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(1,0) = vec( 22 );
         tmp(1,1) = vec( 17 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 15 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );

         try {
            sm = tmp;

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

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) ( 12 ) ( 13 ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( (  7 ) ( 12 ) ( 18 ) ( 14 ) ( 15 ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( -2 ) ( 13 ) ( 22 ) ( 11 ) ( 19 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) ( 15 ) ( 19 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Dense matrix assignment test 6";

      {
         blaze::DynamicMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 15 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 22 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );

         try {
            sm = tmp;

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
         blaze::DynamicMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(2,0) = vec( 22 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 15 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );

         try {
            sm = tmp;

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

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( (  7 ) (    ) (  3 ) (  1 ) ( 12 ) ( 13 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( -2 ) (    ) (  1 ) (  5 ) ( 18 ) ( 14 ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) ( 12 ) ( 18 ) ( 14 ) ( 11 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( 13 ) ( 14 ) ( 22 ) ( 19 ) )
   {
      test_ = "Dense matrix assignment test 7";

      {
         blaze::DynamicMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 11 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 22 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );

         try {
            sm = tmp;

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
         blaze::DynamicMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 22 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );

         try {
            sm = tmp;

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

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) ( 12 ) ( 13 ) ( 19 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) ( 18 ) ( 14 ) ( 11 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( ( 12 ) ( 18 ) ( 14 ) ( 11 ) ( 12 ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( 13 ) ( 14 ) ( 22 ) ( 19 ) ( 14 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( ( 19 ) ( 11 ) ( 12 ) ( 14 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Dense matrix assignment test 8";

      {
         blaze::DynamicMatrix<VT,blaze::rowMajor> tmp( 3UL, 4UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 11 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 22 );
         tmp(1,3) = vec( 19 );
         tmp(2,0) = vec( 19 );
         tmp(2,1) = vec( 11 );
         tmp(2,2) = vec( 12 );
         tmp(2,3) = vec( 14 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );

         try {
            sm = tmp;

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
         blaze::DynamicMatrix<VT,blaze::rowMajor> tmp( 4UL, 3UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(0,2) = vec( 19 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 11 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(2,2) = vec( 12 );
         tmp(3,0) = vec( 22 );
         tmp(3,1) = vec( 19 );
         tmp(3,2) = vec( 14 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );

         try {
            sm = tmp;

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

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( ( 12 ) ( 18 ) ( 14 ) ( 15 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( 18 ) ( 17 ) ( 11 ) ( 19 ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( ( 14 ) ( 11 ) (  3 ) (  1 ) (    ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( 15 ) ( 19 ) (  1 ) (  5 ) (  7 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Sparse matrix assignment test 1";

      {
         blaze::CompressedMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL, 8UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 15 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 17 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 18 ) || sm(0,2) != vec( 14 ) || sm(0,3) != vec( 15 ) ||
             sm(1,0) != vec( 18 ) || sm(1,1) != vec( 17 ) || sm(1,2) != vec( 11 ) || sm(1,3) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 15 ) )\n"
                                        "( ( 18 ) ( 17 ) ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec( 12 )  || sym(0,1) != vec( 18 ) || sym(0,2) != vec( 14 )  || sym(0,3) != vec( 15 )  || sym(0,4) != vec(  5 )  || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( 18 )  || sym(1,1) != vec( 17 ) || sym(1,2) != vec( 11 )  || sym(1,3) != vec( 19 )  || sym(1,4) != vec( -1 )  || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec( 14 )  || sym(2,1) != vec( 11 ) || sym(2,2) != vec(  3 )  || sym(2,3) != vec(  1 )  || !isDefault( sym(2,4) ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( 15 )  || sym(3,1) != vec( 19 ) || sym(3,2) != vec(  1 )  || sym(3,3) != vec(  5 )  || sym(3,4) != vec(  7 )  || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 ) || !isDefault( sym(4,2) ) || sym(4,3) != vec(  7 )  || sym(4,4) != vec(  1 )  || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 )  || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 )  || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 15 ) (  5 ) (    ) )\n"
                                        "( ( 18 ) ( 17 ) ( 11 ) ( 19 ) ( -1 ) (  8 ) )\n"
                                        "( ( 14 ) ( 11 ) (  3 ) (  1 ) (    ) ( -2 ) )\n"
                                        "( ( 15 ) ( 19 ) (  1 ) (  5 ) (  7 ) (    ) )\n"
                                        "( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL, 8UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 17 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 15 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 18 ) ||
             sm(1,0) != vec( 18 ) || sm(1,1) != vec( 17 ) ||
             sm(2,0) != vec( 14 ) || sm(2,1) != vec( 11 ) ||
             sm(3,0) != vec( 15 ) || sm(3,1) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) )\n"
                                        "( ( 18 ) ( 17 ) )\n"
                                        "( ( 14 ) ( 11 ) )\n"
                                        "( ( 15 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec( 12 )  || sym(0,1) != vec( 18 ) || sym(0,2) != vec( 14 )  || sym(0,3) != vec( 15 )  || sym(0,4) != vec(  5 )  || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( 18 )  || sym(1,1) != vec( 17 ) || sym(1,2) != vec( 11 )  || sym(1,3) != vec( 19 )  || sym(1,4) != vec( -1 )  || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec( 14 )  || sym(2,1) != vec( 11 ) || sym(2,2) != vec(  3 )  || sym(2,3) != vec(  1 )  || !isDefault( sym(2,4) ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( 15 )  || sym(3,1) != vec( 19 ) || sym(3,2) != vec(  1 )  || sym(3,3) != vec(  5 )  || sym(3,4) != vec(  7 )  || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 ) || !isDefault( sym(4,2) ) || sym(4,3) != vec(  7 )  || sym(4,4) != vec(  1 )  || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 )  || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 )  || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 15 ) (  5 ) (    ) )\n"
                                        "( ( 18 ) ( 17 ) ( 11 ) ( 19 ) ( -1 ) (  8 ) )\n"
                                        "( ( 14 ) ( 11 ) (  3 ) (  1 ) (    ) ( -2 ) )\n"
                                        "( ( 15 ) ( 19 ) (  1 ) (  5 ) (  7 ) (    ) )\n"
                                        "( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) ( 12 ) ( 13 ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( (  7 ) ( 12 ) ( 18 ) ( 14 ) ( 15 ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( -2 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) ( 15 ) ( 19 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Sparse matrix assignment test 2";

      {
         blaze::CompressedMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL , 8UL);
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 15 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 18 ) || sm(0,2) != vec( 14 ) || sm(0,3) != vec( 15 ) ||
             sm(1,0) != vec( 13 ) || sm(1,1) != vec( 14 ) || sm(1,2) != vec( 11 ) || sm(1,3) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 15 ) )\n"
                                        "( ( 13 ) ( 14 ) ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 ) || sym(0,2) != vec(  7 ) || sym(0,3) != vec( -2 )  || sym(0,4) != vec(  5 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 ) || sym(1,2) != vec( 12 ) || sym(1,3) != vec( 13 )  || sym(1,4) != vec( -1 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec(  7 )  || sym(2,1) != vec( 12 ) || sym(2,2) != vec( 18 ) || sym(2,3) != vec( 14 )  || sym(2,4) != vec( 15 ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( -2 )  || sym(3,1) != vec( 13 ) || sym(3,2) != vec( 14 ) || sym(3,3) != vec( 11 )  || sym(3,4) != vec( 19 ) || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 ) || sym(4,2) != vec( 15 ) || sym(4,3) != vec( 19 )  || sym(4,4) != vec(  1 ) || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 ) || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 ) || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) ( 12 ) ( 13 ) ( -1 ) (  8 ) )\n"
                                        "( (  7 ) ( 12 ) ( 18 ) ( 14 ) ( 15 ) ( -2 ) )\n"
                                        "( ( -2 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) (    ) )\n"
                                        "( (  5 ) ( -1 ) ( 15 ) ( 19 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL, 8UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 15 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 13 ) ||
             sm(1,0) != vec( 18 ) || sm(1,1) != vec( 14 ) ||
             sm(2,0) != vec( 14 ) || sm(2,1) != vec( 11 ) ||
             sm(3,0) != vec( 15 ) || sm(3,1) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 13 ) )\n"
                                        "( ( 18 ) ( 14 ) )\n"
                                        "( ( 14 ) ( 11 ) )\n"
                                        "( ( 15 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 ) || sym(0,2) != vec(  7 ) || sym(0,3) != vec( -2 )  || sym(0,4) != vec(  5 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 ) || sym(1,2) != vec( 12 ) || sym(1,3) != vec( 13 )  || sym(1,4) != vec( -1 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec(  7 )  || sym(2,1) != vec( 12 ) || sym(2,2) != vec( 18 ) || sym(2,3) != vec( 14 )  || sym(2,4) != vec( 15 ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( -2 )  || sym(3,1) != vec( 13 ) || sym(3,2) != vec( 14 ) || sym(3,3) != vec( 11 )  || sym(3,4) != vec( 19 ) || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 ) || sym(4,2) != vec( 15 ) || sym(4,3) != vec( 19 )  || sym(4,4) != vec(  1 ) || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 ) || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 ) || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) ( 12 ) ( 13 ) ( -1 ) (  8 ) )\n"
                                        "( (  7 ) ( 12 ) ( 18 ) ( 14 ) ( 15 ) ( -2 ) )\n"
                                        "( ( -2 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) (    ) )\n"
                                        "( (  5 ) ( -1 ) ( 15 ) ( 19 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( (  7 ) (    ) (  3 ) (  1 ) ( 12 ) ( 13 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( -2 ) (    ) (  1 ) (  5 ) ( 18 ) ( 14 ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) ( 12 ) ( 18 ) ( 14 ) ( 11 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) )
   {
      test_ = "Sparse matrix assignment test 3";

      {
         blaze::CompressedMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL, 8UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 11 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 18 ) || sm(0,2) != vec( 14 ) || sm(0,3) != vec( 11 ) ||
             sm(1,0) != vec( 13 ) || sm(1,1) != vec( 14 ) || sm(1,2) != vec( 11 ) || sm(1,3) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 11 ) )\n"
                                        "( ( 13 ) ( 14 ) ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 )  || sym(0,2) != vec(  7 )  || sym(0,3) != vec( -2 )  || sym(0,4) != vec(  5 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) || !isDefault( sym(1,3) ) || sym(1,4) != vec( -1 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec(  7 )  || !isDefault( sym(2,1) ) || sym(2,2) != vec(  3 )  || sym(2,3) != vec(  1 )  || sym(2,4) != vec( 12 ) || sym(2,5) != vec( 13 )  ||
             sym(3,0) != vec( -2 )  || !isDefault( sym(3,1) ) || sym(3,2) != vec(  1 )  || sym(3,3) != vec(  5 )  || sym(3,4) != vec( 18 ) || sym(3,5) != vec( 14 )  ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 )  || sym(4,2) != vec( 12 )  || sym(4,3) != vec( 18 )  || sym(4,4) != vec( 14 ) || sym(4,5) != vec( 11 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 )  || sym(5,2) != vec( 13 )  || sym(5,3) != vec( 14 )  || sym(5,4) != vec( 11 ) || sym(5,5) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )\n"
                                        "( (  7 ) (    ) (  3 ) (  1 ) ( 12 ) ( 13 ) )\n"
                                        "( ( -2 ) (    ) (  1 ) (  5 ) ( 18 ) ( 14 ) )\n"
                                        "( (  5 ) ( -1 ) ( 12 ) ( 18 ) ( 14 ) ( 11 ) )\n"
                                        "( (    ) (  8 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL, 8UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 11 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 30UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 13 ) ||
             sm(1,0) != vec( 18 ) || sm(1,1) != vec( 14 ) ||
             sm(2,0) != vec( 14 ) || sm(2,1) != vec( 11 ) ||
             sm(3,0) != vec( 11 ) || sm(3,1) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 13 ) )\n"
                                        "( ( 18 ) ( 14 ) )\n"
                                        "( ( 14 ) ( 11 ) )\n"
                                        "( ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 )  || sym(0,2) != vec(  7 )  || sym(0,3) != vec( -2 )  || sym(0,4) != vec(  5 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 )  || !isDefault( sym(1,2) ) || !isDefault( sym(1,3) ) || sym(1,4) != vec( -1 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec(  7 )  || !isDefault( sym(2,1) ) || sym(2,2) != vec(  3 )  || sym(2,3) != vec(  1 )  || sym(2,4) != vec( 12 ) || sym(2,5) != vec( 13 )  ||
             sym(3,0) != vec( -2 )  || !isDefault( sym(3,1) ) || sym(3,2) != vec(  1 )  || sym(3,3) != vec(  5 )  || sym(3,4) != vec( 18 ) || sym(3,5) != vec( 14 )  ||
             sym(4,0) != vec(  5 )  || sym(4,1) != vec( -1 )  || sym(4,2) != vec( 12 )  || sym(4,3) != vec( 18 )  || sym(4,4) != vec( 14 ) || sym(4,5) != vec( 11 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 )  || sym(5,2) != vec( 13 )  || sym(5,3) != vec( 14 )  || sym(5,4) != vec( 11 ) || sym(5,5) != vec( 19 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )\n"
                                        "( (  7 ) (    ) (  3 ) (  1 ) ( 12 ) ( 13 ) )\n"
                                        "( ( -2 ) (    ) (  1 ) (  5 ) ( 18 ) ( 14 ) )\n"
                                        "( (  5 ) ( -1 ) ( 12 ) ( 18 ) ( 14 ) ( 11 ) )\n"
                                        "( (    ) (  8 ) ( 13 ) ( 14 ) ( 11 ) ( 19 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) ( 12 ) ( 13 ) ( 19 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) ( 18 ) ( 14 ) ( 11 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( ( 12 ) ( 18 ) ( 14 ) ( 11 ) ( 12 ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( 13 ) ( 14 ) ( 11 ) ( 19 ) ( 14 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( ( 19 ) ( 11 ) ( 12 ) ( 14 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Sparse matrix assignment test 4";

      {
         blaze::CompressedMatrix<VT,blaze::rowMajor> tmp( 3UL, 4UL, 12UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 11 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );
         tmp(2,0) = vec( 19 );
         tmp(2,1) = vec( 11 );
         tmp(2,2) = vec( 12 );
         tmp(2,3) = vec( 14 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 18 ) || sm(0,2) != vec( 14 ) || sm(0,3) != vec( 11 ) ||
             sm(1,0) != vec( 13 ) || sm(1,1) != vec( 14 ) || sm(1,2) != vec( 11 ) || sm(1,3) != vec( 19 ) ||
             sm(2,0) != vec( 19 ) || sm(2,1) != vec( 11 ) || sm(2,2) != vec( 12 ) || sm(2,3) != vec( 14 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 18 ) ( 14 ) ( 11 ) )\n"
                                        "( ( 13 ) ( 14 ) ( 11 ) ( 19 ) )\n"
                                        "( ( 19 ) ( 11 ) ( 12 ) ( 14 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 ) || sym(0,2) != vec( 12 ) || sym(0,3) != vec( 13 )  || sym(0,4) != vec( 19 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 ) || sym(1,2) != vec( 18 ) || sym(1,3) != vec( 14 )  || sym(1,4) != vec( 11 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec( 12 )  || sym(2,1) != vec( 18 ) || sym(2,2) != vec( 14 ) || sym(2,3) != vec( 11 )  || sym(2,4) != vec( 12 ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( 13 )  || sym(3,1) != vec( 14 ) || sym(3,2) != vec( 11 ) || sym(3,3) != vec( 19 )  || sym(3,4) != vec( 14 ) || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec( 19 )  || sym(4,1) != vec( 11 ) || sym(4,2) != vec( 12 ) || sym(4,3) != vec( 14 )  || sym(4,4) != vec(  1 ) || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 ) || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 ) || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) ( 12 ) ( 13 ) ( 19 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) ( 18 ) ( 14 ) ( 11 ) (  8 ) )\n"
                                        "( ( 12 ) ( 18 ) ( 14 ) ( 11 ) ( 12 ) ( -2 ) )\n"
                                        "( ( 13 ) ( 14 ) ( 11 ) ( 19 ) ( 14 ) (    ) )\n"
                                        "( ( 19 ) ( 11 ) ( 12 ) ( 14 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<VT,blaze::rowMajor> tmp( 4UL, 3UL, 12UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(0,2) = vec( 19 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 11 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(2,2) = vec( 12 );
         tmp(3,0) = vec( 11 );
         tmp(3,1) = vec( 19 );
         tmp(3,2) = vec( 14 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );
         sm = tmp;

         checkRows    ( sym,  6UL );
         checkColumns ( sym,  6UL );
         checkNonZeros( sym, 32UL );

         if( sm(0,0) != vec( 12 ) || sm(0,1) != vec( 13 ) || sm(0,2) != vec( 19 ) ||
             sm(1,0) != vec( 18 ) || sm(1,1) != vec( 14 ) || sm(1,2) != vec( 11 ) ||
             sm(2,0) != vec( 14 ) || sm(2,1) != vec( 11 ) || sm(2,2) != vec( 12 ) ||
             sm(3,0) != vec( 11 ) || sm(3,1) != vec( 19 ) || sm(3,2) != vec( 14 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 12 ) ( 13 ) ( 19 ) )\n"
                                        "( ( 18 ) ( 14 ) ( 11 ) )\n"
                                        "( ( 14 ) ( 11 ) ( 12 ) )\n"
                                        "( ( 11 ) ( 19 ) ( 14 ) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( sym(0,0) != vec(  1 )  || sym(0,1) != vec( -4 ) || sym(0,2) != vec( 12 ) || sym(0,3) != vec( 13 )  || sym(0,4) != vec( 19 ) || !isDefault( sym(0,5) ) ||
             sym(1,0) != vec( -4 )  || sym(1,1) != vec(  2 ) || sym(1,2) != vec( 18 ) || sym(1,3) != vec( 14 )  || sym(1,4) != vec( 11 ) || sym(1,5) != vec(  8 )  ||
             sym(2,0) != vec( 12 )  || sym(2,1) != vec( 18 ) || sym(2,2) != vec( 14 ) || sym(2,3) != vec( 11 )  || sym(2,4) != vec( 12 ) || sym(2,5) != vec( -2 )  ||
             sym(3,0) != vec( 13 )  || sym(3,1) != vec( 14 ) || sym(3,2) != vec( 11 ) || sym(3,3) != vec( 19 )  || sym(3,4) != vec( 14 ) || !isDefault( sym(3,5) ) ||
             sym(4,0) != vec( 19 )  || sym(4,1) != vec( 11 ) || sym(4,2) != vec( 12 ) || sym(4,3) != vec( 14 )  || sym(4,4) != vec(  1 ) || sym(4,5) != vec( -4 )  ||
             !isDefault( sym(5,0) ) || sym(5,1) != vec(  8 ) || sym(5,2) != vec( -2 ) || !isDefault( sym(5,3) ) || sym(5,4) != vec( -4 ) || sym(5,5) != vec(  7 ) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sym << "\n"
                << "   Expected result:\n( (  1 ) ( -4 ) ( 12 ) ( 13 ) ( 19 ) (    ) )\n"
                                        "( ( -4 ) (  2 ) ( 18 ) ( 14 ) ( 11 ) (  8 ) )\n"
                                        "( ( 12 ) ( 18 ) ( 14 ) ( 11 ) ( 12 ) ( -2 ) )\n"
                                        "( ( 13 ) ( 14 ) ( 11 ) ( 19 ) ( 14 ) (    ) )\n"
                                        "( ( 19 ) ( 11 ) ( 12 ) ( 14 ) (  1 ) ( -4 ) )\n"
                                        "( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( ( 12 ) ( 18 ) ( 14 ) ( 15 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( 22 ) ( 17 ) ( 11 ) ( 19 ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( ( 14 ) ( 11 ) (  3 ) (  1 ) (    ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( 15 ) ( 19 ) (  1 ) (  5 ) (  7 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Sparse matrix assignment test 5";

      {
         blaze::CompressedMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL, 8UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 15 );
         tmp(1,0) = vec( 22 );
         tmp(1,1) = vec( 17 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 2UL, 4UL );

         try {
            sm = tmp;

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
         blaze::CompressedMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL, 8UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(1,0) = vec( 22 );
         tmp(1,1) = vec( 17 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 15 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 0UL, 4UL, 2UL );

         try {
            sm = tmp;

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

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) ( 12 ) ( 13 ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( (  7 ) ( 12 ) ( 18 ) ( 14 ) ( 15 ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( -2 ) ( 13 ) ( 22 ) ( 11 ) ( 19 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) ( 15 ) ( 19 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Sparse matrix assignment test 6";

      {
         blaze::CompressedMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL, 8UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 15 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 22 );
         tmp(1,2) = vec( 11 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 1UL, 2UL, 4UL );

         try {
            sm = tmp;

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
         blaze::CompressedMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL, 8UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(2,0) = vec( 22 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 15 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 1UL, 2UL, 4UL, 2UL );

         try {
            sm = tmp;

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

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( (  7 ) (    ) (  3 ) (  1 ) ( 12 ) ( 13 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( -2 ) (    ) (  1 ) (  5 ) ( 18 ) ( 14 ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( (  5 ) ( -1 ) ( 12 ) ( 18 ) ( 14 ) ( 11 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( 13 ) ( 14 ) ( 22 ) ( 19 ) )
   {
      test_ = "Sparse matrix assignment test 7";

      {
         blaze::CompressedMatrix<VT,blaze::rowMajor> tmp( 2UL, 4UL, 8UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 11 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 22 );
         tmp(1,3) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 4UL, 2UL, 2UL, 4UL );

         try {
            sm = tmp;

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
         blaze::CompressedMatrix<VT,blaze::columnMajor> tmp( 4UL, 2UL, 8UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(3,0) = vec( 22 );
         tmp(3,1) = vec( 19 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 4UL, 4UL, 2UL );

         try {
            sm = tmp;

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

   // ( (  1 ) ( -4 ) (  7 ) ( -2 ) (  5 ) (    ) )      ( (  1 ) ( -4 ) ( 12 ) ( 13 ) ( 19 ) (    ) )
   // ( ( -4 ) (  2 ) (    ) (    ) ( -1 ) (  8 ) )      ( ( -4 ) (  2 ) ( 18 ) ( 14 ) ( 11 ) (  8 ) )
   // ( (  7 ) (    ) (  3 ) (  1 ) (    ) ( -2 ) )  =>  ( ( 12 ) ( 18 ) ( 14 ) ( 11 ) ( 12 ) ( -2 ) )
   // ( ( -2 ) (    ) (  1 ) (  5 ) (  7 ) (    ) )      ( ( 13 ) ( 14 ) ( 22 ) ( 19 ) ( 14 ) (    ) )
   // ( (  5 ) ( -1 ) (    ) (  7 ) (  1 ) ( -4 ) )      ( ( 19 ) ( 11 ) ( 12 ) ( 14 ) (  1 ) ( -4 ) )
   // ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )      ( (    ) (  8 ) ( -2 ) (    ) ( -4 ) (  7 ) )
   {
      test_ = "Sparse matrix assignment test 8";

      {
         blaze::CompressedMatrix<VT,blaze::rowMajor> tmp( 3UL, 4UL, 12UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 18 );
         tmp(0,2) = vec( 14 );
         tmp(0,3) = vec( 11 );
         tmp(1,0) = vec( 13 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 22 );
         tmp(1,3) = vec( 19 );
         tmp(2,0) = vec( 19 );
         tmp(2,1) = vec( 11 );
         tmp(2,2) = vec( 12 );
         tmp(2,3) = vec( 14 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 2UL, 0UL, 3UL, 4UL );

         try {
            sm = tmp;

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
         blaze::CompressedMatrix<VT,blaze::rowMajor> tmp( 4UL, 3UL, 12UL );
         tmp(0,0) = vec( 12 );
         tmp(0,1) = vec( 13 );
         tmp(0,2) = vec( 19 );
         tmp(1,0) = vec( 18 );
         tmp(1,1) = vec( 14 );
         tmp(1,2) = vec( 11 );
         tmp(2,0) = vec( 14 );
         tmp(2,1) = vec( 11 );
         tmp(2,2) = vec( 12 );
         tmp(3,0) = vec( 22 );
         tmp(3,1) = vec( 19 );
         tmp(3,2) = vec( 14 );

         ST sym;
         init( sym );

         auto sm = submatrix( sym, 0UL, 2UL, 4UL, 3UL );

         try {
            sm = tmp;

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
void SubmatrixNonScalarTest::checkRows( const Type& matrix, size_t expectedRows ) const
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
void SubmatrixNonScalarTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
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
void SubmatrixNonScalarTest::checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const
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
void SubmatrixNonScalarTest::init( ST& sym )
{
   sym.resize( 6UL );
   sym(0,0) = vec(  1 );
   sym(0,1) = vec( -4 );
   sym(0,2) = vec(  7 );
   sym(0,3) = vec( -2 );
   sym(0,4) = vec(  5 );
   sym(1,1) = vec(  2 );
   sym(1,4) = vec( -1 );
   sym(1,5) = vec(  8 );
   sym(2,2) = vec(  3 );
   sym(2,3) = vec(  1 );
   sym(2,5) = vec( -2 );
   sym(3,3) = vec(  5 );
   sym(3,4) = vec(  7 );
   sym(4,4) = vec(  1 );
   sym(4,5) = vec( -4 );
   sym(5,5) = vec(  7 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setup of a vector.
//
// \param value The value of the vector.
// \return The created vector.
// \exception std::runtime_error Error detected.
//
// This function creates a single vector of size 1. The element of the vector is initialized with
// the given integer value.
*/
inline SubmatrixNonScalarTest::VT SubmatrixNonScalarTest::vec( int value )
{
   return VT( 1UL, value );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the assignment to a submatrix of a non-scalar SymmetricMatrix.
//
// \return void
*/
void runTest()
{
   SubmatrixNonScalarTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the SymmetricMatrix submatrix non-scalar test.
*/
#define RUN_SYMMETRICMATRIX_SUBMATRIXNONSCALAR_TEST \
   blazetest::mathtest::adaptors::symmetricmatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace symmetricmatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
