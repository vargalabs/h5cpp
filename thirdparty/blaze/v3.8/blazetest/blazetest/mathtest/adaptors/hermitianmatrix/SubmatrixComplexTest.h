//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/hermitianmatrix/SubmatrixComplexTest.h
//  \brief Header file for the HermitianMatrix submatrix complex test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_HERMITIANMATRIX_SUBMATRIXCOMPLEXTEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_HERMITIANMATRIX_SUBMATRIXCOMPLEXTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/Submatrix.h>
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
/*!\brief Auxiliary class for assignment tests to a submatrix of a HermitianMatrix.
//
// This class performs assignment tests to a submatrix of a HermitianMatrix with complex element
// type. It performs a series of both compile time as well as runtime tests.
*/
class SubmatrixComplexTest
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
   explicit SubmatrixComplexTest();
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
   template< typename HT > void testAssignment ();
   template< typename HT > void testAddAssign  ();
   template< typename HT > void testSubAssign  ();
   template< typename HT > void testSchurAssign();

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
/*!\brief Test of the assignment to a submatrix of a HermitianMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to a submatrix of a HermitianMatrix. In case an
// error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename HT >  // Type of the Hermitian matrix
void SubmatrixComplexTest::testAssignment()
{
   //=====================================================================================
   // Dense matrix assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 1";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(0,2) = cplx(14,-2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(18, 1);
         mat(1,1) = cplx(17, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) || sm(0,2) != cplx(14,-2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,-2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(1,0) = cplx(18, 1);
         mat(1,1) = cplx(17, 0);
         mat(2,0) = cplx(14, 2);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) ||
             sm(2,0) != cplx(14,2) || sm(2,1) != cplx(11, 1) ||
             sm(3,0) != cplx(15,3) || sm(3,1) != cplx(19, 2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) )\n"
                                        "( (18,1) (17, 0) )\n"
                                        "( (14,2) (11, 1) )\n"
                                        "( (15,3) (19, 2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 2";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 0);
         mat(0,2) = cplx(14, 2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14,-2);
         mat(1,2) = cplx(11, 0);
         mat(1,3) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18, 0) || sm(0,2) != cplx(14,2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,-2) || sm(1,2) != cplx(11,0) || sm(1,3) != cplx(19, 1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18, 0) (14,2) (15,-3) )\n"
                                        "( (13,-2) (14,-2) (11,0) (19, 1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18, 0);
         mat(1,1) = cplx(14, 2);
         mat(2,0) = cplx(14,-2);
         mat(2,1) = cplx(11, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13, 2) ||
             sm(1,0) != cplx(18, 0) || sm(1,1) != cplx(14, 2) ||
             sm(2,0) != cplx(14,-2) || sm(2,1) != cplx(11, 0) ||
             sm(3,0) != cplx(15, 3) || sm(3,1) != cplx(19,-1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13, 2) )\n"
                                        "( (18, 0) (14, 2) )\n"
                                        "( (14,-2) (11, 0) )\n"
                                        "( (15, 3) (19,-1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )
   {
      test_ = "Dense matrix assignment test 3";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 1);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11, 1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18,1) || sm(0,2) != cplx(14, 0) || sm(0,3) != cplx(11,1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18,1) (14, 0) (11,1) )\n"
                                        "( (13,-2) (14,0) (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(11,-1);
         mat(3,1) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13,2) ||
             sm(1,0) != cplx(18,-1) || sm(1,1) != cplx(14,0) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,1) ||
             sm(3,0) != cplx(11,-1) || sm(3,1) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13,2) )\n"
                                        "( (18,-1) (14,0) )\n"
                                        "( (14, 0) (11,1) )\n"
                                        "( (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 4";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11,-1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11, 1);
         mat(1,3) = cplx(19, 0);
         mat(2,0) = cplx(19, 3);
         mat(2,1) = cplx(11, 2);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx(14,-4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(18,-3) || sm(0,2) != cplx(14,0) || sm(0,3) != cplx(11,-1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,1) || sm(1,3) != cplx(19, 0) ||
             sm(2,0) != cplx(19, 3) || sm(2,1) != cplx(11, 2) || sm(2,2) != cplx(12,1) || sm(2,3) != cplx(14,-4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (18,-3) (14, 0) (11,-1) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(13, 2);
         mat(0,2) = cplx(19,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-2);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11,-1);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(11, 1);
         mat(3,1) = cplx(19, 0);
         mat(3,2) = cplx(14, 4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(13, 2) || sm(0,2) != cplx(19,-3) ||
             sm(1,0) != cplx(18, 3) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,-2) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,-1) || sm(2,2) != cplx(12,-1) ||
             sm(3,0) != cplx(11, 1) || sm(3,1) != cplx(19, 0) || sm(3,2) != cplx(14, 4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (13, 2) (19,-3) )\n"
                                        "( (18, 3) (14, 0) (11,-2) )\n"
                                        "( (14, 0) (11,-1) (12,-1) )\n"
                                        "( (11, 1) (19, 0) (14, 4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (22,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 5";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(0,2) = cplx(14,-2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(22, 1);
         mat(1,1) = cplx(17, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(1,0) = cplx(22, 1);
         mat(1,1) = cplx(17, 0);
         mat(2,0) = cplx(14, 2);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (22,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 6";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 0);
         mat(0,2) = cplx(14, 2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(22,-2);
         mat(1,2) = cplx(11, 0);
         mat(1,3) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18, 0);
         mat(1,1) = cplx(14, 2);
         mat(2,0) = cplx(22,-2);
         mat(2,1) = cplx(11, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (22,-1) (19, 0) )
   {
      test_ = "Dense matrix assignment test 7";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 1);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11, 1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(22,-1);
         mat(1,3) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(22,-1);
         mat(3,1) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (22, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 8";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11,-1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(22, 1);
         mat(1,3) = cplx(19, 0);
         mat(2,0) = cplx(19, 3);
         mat(2,1) = cplx(11, 2);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx(14,-4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(13, 2);
         mat(0,2) = cplx(19,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-2);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11,-1);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(22, 1);
         mat(3,1) = cplx(19, 0);
         mat(3,2) = cplx(14, 4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12, 0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,-1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14, 2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15, 3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 9";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(0,2) = cplx(14,-2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(17, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(17, 0);
         mat(2,0) = cplx(14, 2);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12,1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18,0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15,3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 10";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 0);
         mat(0,2) = cplx(14, 2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 2);
         mat(1,2) = cplx(11, 0);
         mat(1,3) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18, 0);
         mat(1,1) = cplx(14, 2);
         mat(2,0) = cplx(14, 2);
         mat(2,1) = cplx(11, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11, 1) (19, 0) )
   {
      test_ = "Dense matrix assignment test 11";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 1);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11, 1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11, 1);
         mat(1,3) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(11, 1);
         mat(3,1) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11,-1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 12";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11,-1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19, 0);
         mat(2,0) = cplx(19, 3);
         mat(2,1) = cplx(11, 2);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx(14,-4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(13, 2);
         mat(0,2) = cplx(19,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-2);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11,-1);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(11,-1);
         mat(3,1) = cplx(19, 0);
         mat(3,2) = cplx(14, 4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 1) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 13";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(0,2) = cplx(14,-2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(18, 1);
         mat(1,1) = cplx(17, 1);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(1,0) = cplx(18, 1);
         mat(1,1) = cplx(17, 1);
         mat(2,0) = cplx(14, 2);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 1) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 14";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 0);
         mat(0,2) = cplx(14, 2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14,-2);
         mat(1,2) = cplx(11, 1);
         mat(1,3) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18, 0);
         mat(1,1) = cplx(14, 2);
         mat(2,0) = cplx(14,-2);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 1) )
   {
      test_ = "Dense matrix assignment test 15";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 1);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11,-1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11, 1);
         mat(1,3) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(11,-1);
         mat(3,1) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 1) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix assignment test 16";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11,-1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11, 1);
         mat(1,3) = cplx(19, 1);
         mat(2,0) = cplx(19, 3);
         mat(2,1) = cplx(11, 2);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx(14,-4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(13, 2);
         mat(0,2) = cplx(19,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-2);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11,-1);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(11, 1);
         mat(3,1) = cplx(19, 1);
         mat(3,2) = cplx(14, 4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Sparse matrix assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 1";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(0,2) = cplx(14,-2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(18, 1);
         mat(1,1) = cplx(17, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) || sm(0,2) != cplx(14,-2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,-2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(1,0) = cplx(18, 1);
         mat(1,1) = cplx(17, 0);
         mat(2,0) = cplx(14, 2);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) ||
             sm(2,0) != cplx(14,2) || sm(2,1) != cplx(11, 1) ||
             sm(3,0) != cplx(15,3) || sm(3,1) != cplx(19, 2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) )\n"
                                        "( (18,1) (17, 0) )\n"
                                        "( (14,2) (11, 1) )\n"
                                        "( (15,3) (19, 2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 2";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 0);
         mat(0,2) = cplx(14, 2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14,-2);
         mat(1,2) = cplx(11, 0);
         mat(1,3) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18, 0) || sm(0,2) != cplx(14,2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,-2) || sm(1,2) != cplx(11,0) || sm(1,3) != cplx(19, 1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18, 0) (14,2) (15,-3) )\n"
                                        "( (13,-2) (14,-2) (11,0) (19, 1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18, 0);
         mat(1,1) = cplx(14, 2);
         mat(2,0) = cplx(14,-2);
         mat(2,1) = cplx(11, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13, 2) ||
             sm(1,0) != cplx(18, 0) || sm(1,1) != cplx(14, 2) ||
             sm(2,0) != cplx(14,-2) || sm(2,1) != cplx(11, 0) ||
             sm(3,0) != cplx(15, 3) || sm(3,1) != cplx(19,-1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13, 2) )\n"
                                        "( (18, 0) (14, 2) )\n"
                                        "( (14,-2) (11, 0) )\n"
                                        "( (15, 3) (19,-1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )
   {
      test_ = "Sparse matrix assignment test 3";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 1);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11, 1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18,1) || sm(0,2) != cplx(14, 0) || sm(0,3) != cplx(11,1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18,1) (14, 0) (11,1) )\n"
                                        "( (13,-2) (14,0) (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(11,-1);
         mat(3,1) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13,2) ||
             sm(1,0) != cplx(18,-1) || sm(1,1) != cplx(14,0) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,1) ||
             sm(3,0) != cplx(11,-1) || sm(3,1) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13,2) )\n"
                                        "( (18,-1) (14,0) )\n"
                                        "( (14, 0) (11,1) )\n"
                                        "( (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 4";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11,-1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11, 1);
         mat(1,3) = cplx(19, 0);
         mat(2,0) = cplx(19, 3);
         mat(2,1) = cplx(11, 2);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx(14,-4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(18,-3) || sm(0,2) != cplx(14,0) || sm(0,3) != cplx(11,-1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,1) || sm(1,3) != cplx(19, 0) ||
             sm(2,0) != cplx(19, 3) || sm(2,1) != cplx(11, 2) || sm(2,2) != cplx(12,1) || sm(2,3) != cplx(14,-4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (18,-3) (14, 0) (11,-1) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(13, 2);
         mat(0,2) = cplx(19,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-2);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11,-1);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(11, 1);
         mat(3,1) = cplx(19, 0);
         mat(3,2) = cplx(14, 4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );
         sm = mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(13, 2) || sm(0,2) != cplx(19,-3) ||
             sm(1,0) != cplx(18, 3) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,-2) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,-1) || sm(2,2) != cplx(12,-1) ||
             sm(3,0) != cplx(11, 1) || sm(3,1) != cplx(19, 0) || sm(3,2) != cplx(14, 4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (13, 2) (19,-3) )\n"
                                        "( (18, 3) (14, 0) (11,-2) )\n"
                                        "( (14, 0) (11,-1) (12,-1) )\n"
                                        "( (11, 1) (19, 0) (14, 4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (22,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 5";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(0,2) = cplx(14,-2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(22, 1);
         mat(1,1) = cplx(17, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(1,0) = cplx(22, 1);
         mat(1,1) = cplx(17, 0);
         mat(2,0) = cplx(14, 2);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (22,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 6";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 0);
         mat(0,2) = cplx(14, 2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(22,-2);
         mat(1,2) = cplx(11, 0);
         mat(1,3) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18, 0);
         mat(1,1) = cplx(14, 2);
         mat(2,0) = cplx(22,-2);
         mat(2,1) = cplx(11, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (22,-1) (19, 0) )
   {
      test_ = "Sparse matrix assignment test 7";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 1);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11, 1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(22,-1);
         mat(1,3) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(22,-1);
         mat(3,1) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (22, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 8";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11,-1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(22, 1);
         mat(1,3) = cplx(19, 0);
         mat(2,0) = cplx(19, 3);
         mat(2,1) = cplx(11, 2);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx(14,-4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(13, 2);
         mat(0,2) = cplx(19,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-2);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11,-1);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(22, 1);
         mat(3,1) = cplx(19, 0);
         mat(3,2) = cplx(14, 4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12, 0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,-1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14, 2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15, 3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 9";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(0,2) = cplx(14,-2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(17, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(17, 0);
         mat(2,0) = cplx(14, 2);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12,1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18,0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15,3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 10";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 0);
         mat(0,2) = cplx(14, 2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 2);
         mat(1,2) = cplx(11, 0);
         mat(1,3) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18, 0);
         mat(1,1) = cplx(14, 2);
         mat(2,0) = cplx(14, 2);
         mat(2,1) = cplx(11, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11, 1) (19, 0) )
   {
      test_ = "Sparse matrix assignment test 11";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 1);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11, 1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11, 1);
         mat(1,3) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(11, 1);
         mat(3,1) = cplx(19, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11,-1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 12";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11,-1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19, 0);
         mat(2,0) = cplx(19, 3);
         mat(2,1) = cplx(11, 2);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx(14,-4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(13, 2);
         mat(0,2) = cplx(19,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-2);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11,-1);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(11,-1);
         mat(3,1) = cplx(19, 0);
         mat(3,2) = cplx(14, 4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 1) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 13";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(0,2) = cplx(14,-2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(18, 1);
         mat(1,1) = cplx(17, 1);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 0);
         mat(0,1) = cplx(18,-1);
         mat(1,0) = cplx(18, 1);
         mat(1,1) = cplx(17, 1);
         mat(2,0) = cplx(14, 2);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 1) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 14";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 0);
         mat(0,2) = cplx(14, 2);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14,-2);
         mat(1,2) = cplx(11, 1);
         mat(1,3) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18, 0);
         mat(1,1) = cplx(14, 2);
         mat(2,0) = cplx(14,-2);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(19,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 1) )
   {
      test_ = "Sparse matrix assignment test 15";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(18, 1);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11,-1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11, 1);
         mat(1,3) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(18,-1);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(11,-1);
         mat(3,1) = cplx(19, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 1) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix assignment test 16";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(11,-1);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11, 1);
         mat(1,3) = cplx(19, 1);
         mat(2,0) = cplx(19, 3);
         mat(2,1) = cplx(11, 2);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx(14,-4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(13, 2);
         mat(0,2) = cplx(19,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(11,-2);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(11,-1);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(11, 1);
         mat(3,1) = cplx(19, 1);
         mat(3,2) = cplx(14, 4);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm = mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the addition assignment to a submatrix of a HermitianMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment to a submatrix of a HermitianMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename HT >  // Type of the Hermitian matrix
void SubmatrixComplexTest::testAddAssign()
{
   //=====================================================================================
   // Dense matrix addition assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 1";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(22, 0);
         mat(0,2) = cplx( 7,-5);
         mat(0,3) = cplx(17,-4);
         mat(1,0) = cplx(22, 0);
         mat(1,1) = cplx(15, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) || sm(0,2) != cplx(14,-2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,-2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(11,0);
         mat(0,1) = cplx(22,0);
         mat(1,0) = cplx(22,0);
         mat(1,1) = cplx(15,0);
         mat(2,0) = cplx( 7,5);
         mat(2,1) = cplx(11,1);
         mat(3,0) = cplx(17,4);
         mat(3,1) = cplx(19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) ||
             sm(2,0) != cplx(14,2) || sm(2,1) != cplx(11, 1) ||
             sm(3,0) != cplx(15,3) || sm(3,1) != cplx(19, 2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) )\n"
                                        "( (18,1) (17, 0) )\n"
                                        "( (14,2) (11, 1) )\n"
                                        "( (15,3) (19, 2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 2";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(15, 0);
         mat(0,2) = cplx(13, 1);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(13,-1);
         mat(1,2) = cplx( 6, 0);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18, 0) || sm(0,2) != cplx(14,2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,-2) || sm(1,2) != cplx(11,0) || sm(1,3) != cplx(19, 1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18, 0) (14,2) (15,-3) )\n"
                                        "( (13,-2) (14,-2) (11,0) (19, 1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(15, 0);
         mat(1,1) = cplx(13, 1);
         mat(2,0) = cplx(13,-1);
         mat(2,1) = cplx( 6, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13, 2) ||
             sm(1,0) != cplx(18, 0) || sm(1,1) != cplx(14, 2) ||
             sm(2,0) != cplx(14,-2) || sm(2,1) != cplx(11, 0) ||
             sm(3,0) != cplx(15, 3) || sm(3,1) != cplx(19,-1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13, 2) )\n"
                                        "( (18, 0) (14, 2) )\n"
                                        "( (14,-2) (11, 0) )\n"
                                        "( (15, 3) (19,-1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )
   {
      test_ = "Dense matrix addition assignment test 3";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(11, 2);
         mat(0,2) = cplx(13, 0);
         mat(0,3) = cplx(15, 1);
         mat(1,0) = cplx(15,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(15,-1);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18,1) || sm(0,2) != cplx(14, 0) || sm(0,3) != cplx(11,1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18,1) (14, 0) (11,1) )\n"
                                        "( (13,-2) (14,0) (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(15, 2);
         mat(1,0) = cplx(11,-2);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(13, 0);
         mat(2,1) = cplx(15, 1);
         mat(3,0) = cplx(15,-1);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13,2) ||
             sm(1,0) != cplx(18,-1) || sm(1,1) != cplx(14,0) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,1) ||
             sm(3,0) != cplx(11,-1) || sm(3,1) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13,2) )\n"
                                        "( (18,-1) (14,0) )\n"
                                        "( (14, 0) (11,1) )\n"
                                        "( (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 4";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx( 5, 4);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(11, 0);
         mat(0,3) = cplx(10,-2);
         mat(1,0) = cplx(15,-1);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(10, 2);
         mat(1,3) = cplx(14, 0);
         mat(2,0) = cplx(14, 3);
         mat(2,1) = cplx(12, 4);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx( 7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(18,-3) || sm(0,2) != cplx(14,0) || sm(0,3) != cplx(11,-1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,1) || sm(1,3) != cplx(19, 0) ||
             sm(2,0) != cplx(19, 3) || sm(2,1) != cplx(11, 2) || sm(2,2) != cplx(12,1) || sm(2,3) != cplx(14,-4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (18,-3) (14, 0) (11,-1) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx( 5,-4);
         mat(0,1) = cplx(15, 1);
         mat(0,2) = cplx(14,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(12,-4);
         mat(2,0) = cplx(11, 0);
         mat(2,1) = cplx(10,-2);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(10, 2);
         mat(3,1) = cplx(14, 0);
         mat(3,2) = cplx( 7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(13, 2) || sm(0,2) != cplx(19,-3) ||
             sm(1,0) != cplx(18, 3) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,-2) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,-1) || sm(2,2) != cplx(12,-1) ||
             sm(3,0) != cplx(11, 1) || sm(3,1) != cplx(19, 0) || sm(3,2) != cplx(14, 4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (13, 2) (19,-3) )\n"
                                        "( (18, 3) (14, 0) (11,-2) )\n"
                                        "( (14, 0) (11,-1) (12,-1) )\n"
                                        "( (11, 1) (19, 0) (14, 4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (22,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 5";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(22, 0);
         mat(0,2) = cplx( 7,-5);
         mat(0,3) = cplx(17,-4);
         mat(1,0) = cplx(26, 0);
         mat(1,1) = cplx(15, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(11,0);
         mat(0,1) = cplx(22,0);
         mat(1,0) = cplx(26,0);
         mat(1,1) = cplx(15,0);
         mat(2,0) = cplx( 7,5);
         mat(2,1) = cplx(11,1);
         mat(3,0) = cplx(17,4);
         mat(3,1) = cplx(19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (22,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 6";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(15, 0);
         mat(0,2) = cplx(13, 1);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(21,-1);
         mat(1,2) = cplx( 6, 0);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(15, 0);
         mat(1,1) = cplx(13, 1);
         mat(2,0) = cplx(21,-1);
         mat(2,1) = cplx( 6, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (22,-1) (19, 0) )
   {
      test_ = "Dense matrix addition assignment test 7";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(11, 2);
         mat(0,2) = cplx(13, 0);
         mat(0,3) = cplx(15, 1);
         mat(1,0) = cplx(15,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(26,-1);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(15, 2);
         mat(1,0) = cplx(11,-2);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(13, 0);
         mat(2,1) = cplx(15, 1);
         mat(3,0) = cplx(26,-1);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (22, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 8";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx( 5, 4);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(11, 0);
         mat(0,3) = cplx(10,-2);
         mat(1,0) = cplx(15,-1);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(21, 2);
         mat(1,3) = cplx(14, 0);
         mat(2,0) = cplx(14, 3);
         mat(2,1) = cplx(12, 4);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx( 7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx( 5,-4);
         mat(0,1) = cplx(15, 1);
         mat(0,2) = cplx(14,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(12,-4);
         mat(2,0) = cplx(11, 0);
         mat(2,1) = cplx(10,-2);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(21, 2);
         mat(3,1) = cplx(14, 0);
         mat(3,2) = cplx( 7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12, 0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,-1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14, 2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15, 3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 9";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(22, 0);
         mat(0,2) = cplx( 7,-5);
         mat(0,3) = cplx(17,-4);
         mat(1,0) = cplx(22,-2);
         mat(1,1) = cplx(15, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(22, 0);
         mat(1,0) = cplx(22,-2);
         mat(1,1) = cplx(15, 0);
         mat(2,0) = cplx( 7, 5);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(17, 4);
         mat(3,1) = cplx(19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12,1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18,0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15,3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 10";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(15, 0);
         mat(0,2) = cplx(13, 1);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(13,-3);
         mat(1,2) = cplx( 6, 0);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(15, 0);
         mat(1,1) = cplx(13, 1);
         mat(2,0) = cplx(13,-3);
         mat(2,1) = cplx( 6, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11, 1) (19, 0) )
   {
      test_ = "Dense matrix addition assignment test 11";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(11, 2);
         mat(0,2) = cplx(13, 0);
         mat(0,3) = cplx(15, 1);
         mat(1,0) = cplx(15,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(15, 1);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(15, 2);
         mat(1,0) = cplx(11,-2);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(13, 0);
         mat(2,1) = cplx(15, 1);
         mat(3,0) = cplx(15, 1);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11,-1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 12";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx( 5, 4);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(11, 0);
         mat(0,3) = cplx(10,-2);
         mat(1,0) = cplx(15,-1);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(10, 0);
         mat(1,3) = cplx(14, 0);
         mat(2,0) = cplx(14, 3);
         mat(2,1) = cplx(12, 4);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx( 7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx( 5,-4);
         mat(0,1) = cplx(15, 1);
         mat(0,2) = cplx(14,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(12,-4);
         mat(2,0) = cplx(11, 0);
         mat(2,1) = cplx(10,-2);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(10, 0);
         mat(3,1) = cplx(14, 0);
         mat(3,2) = cplx( 7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 1) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 13";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(22, 0);
         mat(0,2) = cplx( 7,-5);
         mat(0,3) = cplx(17,-4);
         mat(1,0) = cplx(22, 0);
         mat(1,1) = cplx(15, 1);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(11,0);
         mat(0,1) = cplx(22,0);
         mat(1,0) = cplx(22,0);
         mat(1,1) = cplx(15,1);
         mat(2,0) = cplx( 7,5);
         mat(2,1) = cplx(11,1);
         mat(3,0) = cplx(17,4);
         mat(3,1) = cplx(19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 1) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 14";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(15, 0);
         mat(0,2) = cplx(13, 1);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(13,-1);
         mat(1,2) = cplx( 6, 1);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(15, 0);
         mat(1,1) = cplx(13, 1);
         mat(2,0) = cplx(13,-1);
         mat(2,1) = cplx( 6, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 1) )
   {
      test_ = "Dense matrix addition assignment test 15";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(11, 2);
         mat(0,2) = cplx(13, 0);
         mat(0,3) = cplx(15, 1);
         mat(1,0) = cplx(15,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(15,-1);
         mat(1,3) = cplx(12, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(15, 2);
         mat(1,0) = cplx(11,-2);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(13, 0);
         mat(2,1) = cplx(15, 1);
         mat(3,0) = cplx(15,-1);
         mat(3,1) = cplx(12, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 1) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix addition assignment test 16";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx( 5, 4);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(11, 0);
         mat(0,3) = cplx(10,-2);
         mat(1,0) = cplx(15,-1);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(10, 2);
         mat(1,3) = cplx(14, 1);
         mat(2,0) = cplx(14, 3);
         mat(2,1) = cplx(12, 4);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx( 7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx( 5,-4);
         mat(0,1) = cplx(15, 1);
         mat(0,2) = cplx(14,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 1);
         mat(1,2) = cplx(12,-4);
         mat(2,0) = cplx(11, 0);
         mat(2,1) = cplx(10,-2);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(10, 2);
         mat(3,1) = cplx(14, 1);
         mat(3,2) = cplx( 7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Sparse matrix addition assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 1";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(22, 0);
         mat(0,2) = cplx( 7,-5);
         mat(0,3) = cplx(17,-4);
         mat(1,0) = cplx(22, 0);
         mat(1,1) = cplx(15, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) || sm(0,2) != cplx(14,-2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,-2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(11,0);
         mat(0,1) = cplx(22,0);
         mat(1,0) = cplx(22,0);
         mat(1,1) = cplx(15,0);
         mat(2,0) = cplx( 7,5);
         mat(2,1) = cplx(11,1);
         mat(3,0) = cplx(17,4);
         mat(3,1) = cplx(19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) ||
             sm(2,0) != cplx(14,2) || sm(2,1) != cplx(11, 1) ||
             sm(3,0) != cplx(15,3) || sm(3,1) != cplx(19, 2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) )\n"
                                        "( (18,1) (17, 0) )\n"
                                        "( (14,2) (11, 1) )\n"
                                        "( (15,3) (19, 2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 2";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(15, 0);
         mat(0,2) = cplx(13, 1);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(13,-1);
         mat(1,2) = cplx( 6, 0);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18, 0) || sm(0,2) != cplx(14,2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,-2) || sm(1,2) != cplx(11,0) || sm(1,3) != cplx(19, 1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18, 0) (14,2) (15,-3) )\n"
                                        "( (13,-2) (14,-2) (11,0) (19, 1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(15, 0);
         mat(1,1) = cplx(13, 1);
         mat(2,0) = cplx(13,-1);
         mat(2,1) = cplx( 6, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13, 2) ||
             sm(1,0) != cplx(18, 0) || sm(1,1) != cplx(14, 2) ||
             sm(2,0) != cplx(14,-2) || sm(2,1) != cplx(11, 0) ||
             sm(3,0) != cplx(15, 3) || sm(3,1) != cplx(19,-1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13, 2) )\n"
                                        "( (18, 0) (14, 2) )\n"
                                        "( (14,-2) (11, 0) )\n"
                                        "( (15, 3) (19,-1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )
   {
      test_ = "Sparse matrix addition assignment test 3";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(11, 2);
         mat(0,2) = cplx(13, 0);
         mat(0,3) = cplx(15, 1);
         mat(1,0) = cplx(15,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(15,-1);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18,1) || sm(0,2) != cplx(14, 0) || sm(0,3) != cplx(11,1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18,1) (14, 0) (11,1) )\n"
                                        "( (13,-2) (14,0) (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(15, 2);
         mat(1,0) = cplx(11,-2);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(13, 0);
         mat(2,1) = cplx(15, 1);
         mat(3,0) = cplx(15,-1);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13,2) ||
             sm(1,0) != cplx(18,-1) || sm(1,1) != cplx(14,0) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,1) ||
             sm(3,0) != cplx(11,-1) || sm(3,1) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13,2) )\n"
                                        "( (18,-1) (14,0) )\n"
                                        "( (14, 0) (11,1) )\n"
                                        "( (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 4";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx( 5, 4);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(11, 0);
         mat(0,3) = cplx(10,-2);
         mat(1,0) = cplx(15,-1);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(10, 2);
         mat(1,3) = cplx(14, 0);
         mat(2,0) = cplx(14, 3);
         mat(2,1) = cplx(12, 4);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx( 7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(18,-3) || sm(0,2) != cplx(14,0) || sm(0,3) != cplx(11,-1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,1) || sm(1,3) != cplx(19, 0) ||
             sm(2,0) != cplx(19, 3) || sm(2,1) != cplx(11, 2) || sm(2,2) != cplx(12,1) || sm(2,3) != cplx(14,-4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (18,-3) (14, 0) (11,-1) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx( 5,-4);
         mat(0,1) = cplx(15, 1);
         mat(0,2) = cplx(14,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(12,-4);
         mat(2,0) = cplx(11, 0);
         mat(2,1) = cplx(10,-2);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(10, 2);
         mat(3,1) = cplx(14, 0);
         mat(3,2) = cplx( 7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );
         sm += mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(13, 2) || sm(0,2) != cplx(19,-3) ||
             sm(1,0) != cplx(18, 3) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,-2) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,-1) || sm(2,2) != cplx(12,-1) ||
             sm(3,0) != cplx(11, 1) || sm(3,1) != cplx(19, 0) || sm(3,2) != cplx(14, 4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (13, 2) (19,-3) )\n"
                                        "( (18, 3) (14, 0) (11,-2) )\n"
                                        "( (14, 0) (11,-1) (12,-1) )\n"
                                        "( (11, 1) (19, 0) (14, 4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (22,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 5";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(22, 0);
         mat(0,2) = cplx( 7,-5);
         mat(0,3) = cplx(17,-4);
         mat(1,0) = cplx(26, 0);
         mat(1,1) = cplx(15, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(11,0);
         mat(0,1) = cplx(22,0);
         mat(1,0) = cplx(26,0);
         mat(1,1) = cplx(15,0);
         mat(2,0) = cplx( 7,5);
         mat(2,1) = cplx(11,1);
         mat(3,0) = cplx(17,4);
         mat(3,1) = cplx(19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (22,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 6";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(15, 0);
         mat(0,2) = cplx(13, 1);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(21,-1);
         mat(1,2) = cplx( 6, 0);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(15, 0);
         mat(1,1) = cplx(13, 1);
         mat(2,0) = cplx(21,-1);
         mat(2,1) = cplx( 6, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (22,-1) (19, 0) )
   {
      test_ = "Sparse matrix addition assignment test 7";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(11, 2);
         mat(0,2) = cplx(13, 0);
         mat(0,3) = cplx(15, 1);
         mat(1,0) = cplx(15,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(26,-1);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(15, 2);
         mat(1,0) = cplx(11,-2);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(13, 0);
         mat(2,1) = cplx(15, 1);
         mat(3,0) = cplx(26,-1);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (22, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 8";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx( 5, 4);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(11, 0);
         mat(0,3) = cplx(10,-2);
         mat(1,0) = cplx(15,-1);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(21, 2);
         mat(1,3) = cplx(14, 0);
         mat(2,0) = cplx(14, 3);
         mat(2,1) = cplx(12, 4);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx( 7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx( 5,-4);
         mat(0,1) = cplx(15, 1);
         mat(0,2) = cplx(14,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(12,-4);
         mat(2,0) = cplx(11, 0);
         mat(2,1) = cplx(10,-2);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(21, 2);
         mat(3,1) = cplx(14, 0);
         mat(3,2) = cplx( 7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12, 0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,-1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14, 2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15, 3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 9";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 12UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(22, 0);
         mat(0,2) = cplx( 7,-5);
         mat(0,3) = cplx(17,-4);
         mat(1,0) = cplx(22,-2);
         mat(1,1) = cplx(15, 0);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(22, 0);
         mat(1,0) = cplx(22,-2);
         mat(1,1) = cplx(15, 0);
         mat(2,0) = cplx( 7, 5);
         mat(2,1) = cplx(11, 1);
         mat(3,0) = cplx(17, 4);
         mat(3,1) = cplx(19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12,1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18,0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15,3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 10";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(15, 0);
         mat(0,2) = cplx(13, 1);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(13,-3);
         mat(1,2) = cplx( 6, 0);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(15, 0);
         mat(1,1) = cplx(13, 1);
         mat(2,0) = cplx(13,-3);
         mat(2,1) = cplx( 6, 0);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11, 1) (19, 0) )
   {
      test_ = "Sparse matrix addition assignment test 11";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(11, 2);
         mat(0,2) = cplx(13, 0);
         mat(0,3) = cplx(15, 1);
         mat(1,0) = cplx(15,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(15, 1);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(15, 2);
         mat(1,0) = cplx(11,-2);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(13, 0);
         mat(2,1) = cplx(15, 1);
         mat(3,0) = cplx(15, 1);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11,-1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 12";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx( 5, 4);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(11, 0);
         mat(0,3) = cplx(10,-2);
         mat(1,0) = cplx(15,-1);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(10, 0);
         mat(1,3) = cplx(14, 0);
         mat(2,0) = cplx(14, 3);
         mat(2,1) = cplx(12, 4);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx( 7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx( 5,-4);
         mat(0,1) = cplx(15, 1);
         mat(0,2) = cplx(14,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(12,-4);
         mat(2,0) = cplx(11, 0);
         mat(2,1) = cplx(10,-2);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(10, 0);
         mat(3,1) = cplx(14, 0);
         mat(3,2) = cplx( 7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 1) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 13";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(22, 0);
         mat(0,2) = cplx( 7,-5);
         mat(0,3) = cplx(17,-4);
         mat(1,0) = cplx(22, 0);
         mat(1,1) = cplx(15, 1);
         mat(1,2) = cplx(11,-1);
         mat(1,3) = cplx(19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(11,0);
         mat(0,1) = cplx(22,0);
         mat(1,0) = cplx(22,0);
         mat(1,1) = cplx(15,1);
         mat(2,0) = cplx( 7,5);
         mat(2,1) = cplx(11,1);
         mat(3,0) = cplx(17,4);
         mat(3,1) = cplx(19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 1) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 14";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(15, 0);
         mat(0,2) = cplx(13, 1);
         mat(0,3) = cplx(15,-3);
         mat(1,0) = cplx(13,-2);
         mat(1,1) = cplx(13,-1);
         mat(1,2) = cplx( 6, 1);
         mat(1,3) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(13, 2);
         mat(1,0) = cplx(15, 0);
         mat(1,1) = cplx(13, 1);
         mat(2,0) = cplx(13,-1);
         mat(2,1) = cplx( 6, 1);
         mat(3,0) = cplx(15, 3);
         mat(3,1) = cplx(12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 1) )
   {
      test_ = "Sparse matrix addition assignment test 15";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(12,-1);
         mat(0,1) = cplx(11, 2);
         mat(0,2) = cplx(13, 0);
         mat(0,3) = cplx(15, 1);
         mat(1,0) = cplx(15,-2);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(15,-1);
         mat(1,3) = cplx(12, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(12, 1);
         mat(0,1) = cplx(15, 2);
         mat(1,0) = cplx(11,-2);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(13, 0);
         mat(2,1) = cplx(15, 1);
         mat(3,0) = cplx(15,-1);
         mat(3,1) = cplx(12, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 1) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix addition assignment test 16";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx( 5, 4);
         mat(0,1) = cplx(18,-3);
         mat(0,2) = cplx(11, 0);
         mat(0,3) = cplx(10,-2);
         mat(1,0) = cplx(15,-1);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx(10, 2);
         mat(1,3) = cplx(14, 1);
         mat(2,0) = cplx(14, 3);
         mat(2,1) = cplx(12, 4);
         mat(2,2) = cplx(12, 1);
         mat(2,3) = cplx( 7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx( 5,-4);
         mat(0,1) = cplx(15, 1);
         mat(0,2) = cplx(14,-3);
         mat(1,0) = cplx(18, 3);
         mat(1,1) = cplx(14, 1);
         mat(1,2) = cplx(12,-4);
         mat(2,0) = cplx(11, 0);
         mat(2,1) = cplx(10,-2);
         mat(2,2) = cplx(12,-1);
         mat(3,0) = cplx(10, 2);
         mat(3,1) = cplx(14, 1);
         mat(3,2) = cplx( 7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm += mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the subtraction assignment to a submatrix of a HermitianMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment to a submatrix of a HermitianMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename HT >  // Type of the Hermitian matrix
void SubmatrixComplexTest::testSubAssign()
{
   //=====================================================================================
   // Dense matrix subtraction assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 1";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-11,0);
         mat(0,1) = cplx(-22,0);
         mat(0,2) = cplx( -7,5);
         mat(0,3) = cplx(-17,4);
         mat(1,0) = cplx(-22,0);
         mat(1,1) = cplx(-15,0);
         mat(1,2) = cplx(-11,1);
         mat(1,3) = cplx(-19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) || sm(0,2) != cplx(14,-2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,-2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-11, 0);
         mat(0,1) = cplx(-22, 0);
         mat(1,0) = cplx(-22, 0);
         mat(1,1) = cplx(-15, 0);
         mat(2,0) = cplx( -7,-5);
         mat(2,1) = cplx(-11,-1);
         mat(3,0) = cplx(-17,-4);
         mat(3,1) = cplx(-19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) ||
             sm(2,0) != cplx(14,2) || sm(2,1) != cplx(11, 1) ||
             sm(3,0) != cplx(15,3) || sm(3,1) != cplx(19, 2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) )\n"
                                        "( (18,1) (17, 0) )\n"
                                        "( (14,2) (11, 1) )\n"
                                        "( (15,3) (19, 2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 2";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-15, 0);
         mat(0,2) = cplx(-13,-1);
         mat(0,3) = cplx(-15, 3);
         mat(1,0) = cplx(-13, 2);
         mat(1,1) = cplx(-13, 1);
         mat(1,2) = cplx( -6, 0);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18, 0) || sm(0,2) != cplx(14,2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,-2) || sm(1,2) != cplx(11,0) || sm(1,3) != cplx(19, 1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18, 0) (14,2) (15,-3) )\n"
                                        "( (13,-2) (14,-2) (11,0) (19, 1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-13,-2);
         mat(1,0) = cplx(-15, 0);
         mat(1,1) = cplx(-13,-1);
         mat(2,0) = cplx(-13, 1);
         mat(2,1) = cplx( -6, 0);
         mat(3,0) = cplx(-15,-3);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13, 2) ||
             sm(1,0) != cplx(18, 0) || sm(1,1) != cplx(14, 2) ||
             sm(2,0) != cplx(14,-2) || sm(2,1) != cplx(11, 0) ||
             sm(3,0) != cplx(15, 3) || sm(3,1) != cplx(19,-1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13, 2) )\n"
                                        "( (18, 0) (14, 2) )\n"
                                        "( (14,-2) (11, 0) )\n"
                                        "( (15, 3) (19,-1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 3";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-11,-2);
         mat(0,2) = cplx(-13, 0);
         mat(0,3) = cplx(-15,-1);
         mat(1,0) = cplx(-15, 2);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-15, 1);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18,1) || sm(0,2) != cplx(14, 0) || sm(0,3) != cplx(11,1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18,1) (14, 0) (11,1) )\n"
                                        "( (13,-2) (14,0) (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-15,-2);
         mat(1,0) = cplx(-11, 2);
         mat(1,1) = cplx(-14, 0);
         mat(2,0) = cplx(-13, 0);
         mat(2,1) = cplx(-15,-1);
         mat(3,0) = cplx(-15, 1);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13,2) ||
             sm(1,0) != cplx(18,-1) || sm(1,1) != cplx(14,0) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,1) ||
             sm(3,0) != cplx(11,-1) || sm(3,1) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13,2) )\n"
                                        "( (18,-1) (14,0) )\n"
                                        "( (14, 0) (11,1) )\n"
                                        "( (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 4";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx( -5,-4);
         mat(0,1) = cplx(-18, 3);
         mat(0,2) = cplx(-11, 0);
         mat(0,3) = cplx(-10, 2);
         mat(1,0) = cplx(-15, 1);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-10,-2);
         mat(1,3) = cplx(-14, 0);
         mat(2,0) = cplx(-14,-3);
         mat(2,1) = cplx(-12,-4);
         mat(2,2) = cplx(-12,-1);
         mat(2,3) = cplx( -7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(18,-3) || sm(0,2) != cplx(14,0) || sm(0,3) != cplx(11,-1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,1) || sm(1,3) != cplx(19, 0) ||
             sm(2,0) != cplx(19, 3) || sm(2,1) != cplx(11, 2) || sm(2,2) != cplx(12,1) || sm(2,3) != cplx(14,-4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (18,-3) (14, 0) (11,-1) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx( -5, 4);
         mat(0,1) = cplx(-15,-1);
         mat(0,2) = cplx(-14, 3);
         mat(1,0) = cplx(-18,-3);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-12, 4);
         mat(2,0) = cplx(-11, 0);
         mat(2,1) = cplx(-10, 2);
         mat(2,2) = cplx(-12, 1);
         mat(3,0) = cplx(-10,-2);
         mat(3,1) = cplx(-14, 0);
         mat(3,2) = cplx( -7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(13, 2) || sm(0,2) != cplx(19,-3) ||
             sm(1,0) != cplx(18, 3) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,-2) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,-1) || sm(2,2) != cplx(12,-1) ||
             sm(3,0) != cplx(11, 1) || sm(3,1) != cplx(19, 0) || sm(3,2) != cplx(14, 4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (13, 2) (19,-3) )\n"
                                        "( (18, 3) (14, 0) (11,-2) )\n"
                                        "( (14, 0) (11,-1) (12,-1) )\n"
                                        "( (11, 1) (19, 0) (14, 4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (22,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 5";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-11,0);
         mat(0,1) = cplx(-22,0);
         mat(0,2) = cplx( -7,5);
         mat(0,3) = cplx(-17,4);
         mat(1,0) = cplx(-26,0);
         mat(1,1) = cplx(-15,0);
         mat(1,2) = cplx(-11,1);
         mat(1,3) = cplx(-19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-11, 0);
         mat(0,1) = cplx(-22, 0);
         mat(1,0) = cplx(-26, 0);
         mat(1,1) = cplx(-15, 0);
         mat(2,0) = cplx( -7,-5);
         mat(2,1) = cplx(-11,-1);
         mat(3,0) = cplx(-17,-4);
         mat(3,1) = cplx(-19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (22,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 6";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-15, 0);
         mat(0,2) = cplx(-13,-1);
         mat(0,3) = cplx(-15, 3);
         mat(1,0) = cplx(-13, 2);
         mat(1,1) = cplx(-21, 1);
         mat(1,2) = cplx( -6, 0);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-13,-2);
         mat(1,0) = cplx(-15, 0);
         mat(1,1) = cplx(-13,-1);
         mat(2,0) = cplx(-21, 1);
         mat(2,1) = cplx( -6, 0);
         mat(3,0) = cplx(-15,-3);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (22,-1) (19, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 7";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-11,-2);
         mat(0,2) = cplx(-13, 0);
         mat(0,3) = cplx(-15,-1);
         mat(1,0) = cplx(-15, 2);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-26, 1);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-15,-2);
         mat(1,0) = cplx(-11, 2);
         mat(1,1) = cplx(-14, 0);
         mat(2,0) = cplx(-13, 0);
         mat(2,1) = cplx(-15,-1);
         mat(3,0) = cplx(-26, 1);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (22, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 8";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx( -5,-4);
         mat(0,1) = cplx(-18, 3);
         mat(0,2) = cplx(-11, 0);
         mat(0,3) = cplx(-10, 2);
         mat(1,0) = cplx(-15, 1);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-21,-2);
         mat(1,3) = cplx(-14, 0);
         mat(2,0) = cplx(-14,-3);
         mat(2,1) = cplx(-12,-4);
         mat(2,2) = cplx(-12,-1);
         mat(2,3) = cplx( -7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx( -5, 4);
         mat(0,1) = cplx(-15,-1);
         mat(0,2) = cplx(-14, 3);
         mat(1,0) = cplx(-18,-3);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-12, 4);
         mat(2,0) = cplx(-11, 0);
         mat(2,1) = cplx(-10, 2);
         mat(2,2) = cplx(-12, 1);
         mat(3,0) = cplx(-21,-2);
         mat(3,1) = cplx(-14, 0);
         mat(3,2) = cplx( -7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12, 0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,-1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14, 2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15, 3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 9";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-11,0);
         mat(0,1) = cplx(-22,0);
         mat(0,2) = cplx( -7,5);
         mat(0,3) = cplx(-17,4);
         mat(1,0) = cplx(-22,2);
         mat(1,1) = cplx(-15,0);
         mat(1,2) = cplx(-11,1);
         mat(1,3) = cplx(-19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-11, 0);
         mat(0,1) = cplx(-22, 0);
         mat(1,0) = cplx(-22, 2);
         mat(1,1) = cplx(-15, 0);
         mat(2,0) = cplx( -7,-5);
         mat(2,1) = cplx(-11,-1);
         mat(3,0) = cplx(-17,-4);
         mat(3,1) = cplx(-19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12,1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18,0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15,3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 10";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-15, 0);
         mat(0,2) = cplx(-13,-1);
         mat(0,3) = cplx(-15, 3);
         mat(1,0) = cplx(-13, 2);
         mat(1,1) = cplx(-13, 3);
         mat(1,2) = cplx( -6, 0);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-13,-2);
         mat(1,0) = cplx(-15, 0);
         mat(1,1) = cplx(-13,-1);
         mat(2,0) = cplx(-13, 3);
         mat(2,1) = cplx( -6, 0);
         mat(3,0) = cplx(-15,-3);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11, 1) (19, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 11";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-11,-2);
         mat(0,2) = cplx(-13, 0);
         mat(0,3) = cplx(-15,-1);
         mat(1,0) = cplx(-15, 2);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-15,-1);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-15,-2);
         mat(1,0) = cplx(-11, 2);
         mat(1,1) = cplx(-14, 0);
         mat(2,0) = cplx(-13, 0);
         mat(2,1) = cplx(-15,-1);
         mat(3,0) = cplx(-15,-1);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11,-1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 12";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx( -5,-4);
         mat(0,1) = cplx(-18, 3);
         mat(0,2) = cplx(-11, 0);
         mat(0,3) = cplx(-10, 2);
         mat(1,0) = cplx(-15, 1);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-10, 0);
         mat(1,3) = cplx(-14, 0);
         mat(2,0) = cplx(-14,-3);
         mat(2,1) = cplx(-12,-4);
         mat(2,2) = cplx(-12,-1);
         mat(2,3) = cplx( -7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx( -5, 4);
         mat(0,1) = cplx(-15,-1);
         mat(0,2) = cplx(-14, 3);
         mat(1,0) = cplx(-18,-3);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-12, 4);
         mat(2,0) = cplx(-11, 0);
         mat(2,1) = cplx(-10, 2);
         mat(2,2) = cplx(-12, 1);
         mat(3,0) = cplx(-10, 0);
         mat(3,1) = cplx(-14, 0);
         mat(3,2) = cplx( -7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 1) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 13";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-11, 0);
         mat(0,1) = cplx(-22, 0);
         mat(0,2) = cplx( -7, 5);
         mat(0,3) = cplx(-17, 4);
         mat(1,0) = cplx(-22, 0);
         mat(1,1) = cplx(-15,-1);
         mat(1,2) = cplx(-11, 1);
         mat(1,3) = cplx(-19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-11, 0);
         mat(0,1) = cplx(-22, 0);
         mat(1,0) = cplx(-22, 0);
         mat(1,1) = cplx(-15,-1);
         mat(2,0) = cplx( -7,-5);
         mat(2,1) = cplx(-11,-1);
         mat(3,0) = cplx(-17,-4);
         mat(3,1) = cplx(-19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 1) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 14";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-15, 0);
         mat(0,2) = cplx(-13,-1);
         mat(0,3) = cplx(-15, 3);
         mat(1,0) = cplx(-13, 2);
         mat(1,1) = cplx(-13, 1);
         mat(1,2) = cplx( -6,-1);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-13,-2);
         mat(1,0) = cplx(-15, 0);
         mat(1,1) = cplx(-13,-1);
         mat(2,0) = cplx(-13, 1);
         mat(2,1) = cplx( -6,-1);
         mat(3,0) = cplx(-15,-3);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 1) )
   {
      test_ = "Dense matrix subtraction assignment test 15";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-11,-2);
         mat(0,2) = cplx(-13, 0);
         mat(0,3) = cplx(-15,-1);
         mat(1,0) = cplx(-15, 2);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-15, 1);
         mat(1,3) = cplx(-12,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-15,-2);
         mat(1,0) = cplx(-11, 2);
         mat(1,1) = cplx(-14, 0);
         mat(2,0) = cplx(-13, 0);
         mat(2,1) = cplx(-15,-1);
         mat(3,0) = cplx(-15, 1);
         mat(3,1) = cplx(-12,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 1) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix subtraction assignment test 16";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx( -5,-4);
         mat(0,1) = cplx(-18, 3);
         mat(0,2) = cplx(-11, 0);
         mat(0,3) = cplx(-10, 2);
         mat(1,0) = cplx(-15, 1);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-10,-2);
         mat(1,3) = cplx(-14,-1);
         mat(2,0) = cplx(-14,-3);
         mat(2,1) = cplx(-12,-4);
         mat(2,2) = cplx(-12,-1);
         mat(2,3) = cplx( -7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx( -5, 4);
         mat(0,1) = cplx(-15,-1);
         mat(0,2) = cplx(-14, 3);
         mat(1,0) = cplx(-18,-3);
         mat(1,1) = cplx(-14,-1);
         mat(1,2) = cplx(-12, 4);
         mat(2,0) = cplx(-11, 0);
         mat(2,1) = cplx(-10, 2);
         mat(2,2) = cplx(-12, 1);
         mat(3,0) = cplx(-10,-2);
         mat(3,1) = cplx(-14,-1);
         mat(3,2) = cplx( -7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Sparse matrix subtraction assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 1";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-11,0);
         mat(0,1) = cplx(-22,0);
         mat(0,2) = cplx( -7,5);
         mat(0,3) = cplx(-17,4);
         mat(1,0) = cplx(-22,0);
         mat(1,1) = cplx(-15,0);
         mat(1,2) = cplx(-11,1);
         mat(1,3) = cplx(-19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) || sm(0,2) != cplx(14,-2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,-2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-11, 0);
         mat(0,1) = cplx(-22, 0);
         mat(1,0) = cplx(-22, 0);
         mat(1,1) = cplx(-15, 0);
         mat(2,0) = cplx( -7,-5);
         mat(2,1) = cplx(-11,-1);
         mat(3,0) = cplx(-17,-4);
         mat(3,1) = cplx(-19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,0) || sm(0,1) != cplx(18,-1) ||
             sm(1,0) != cplx(18,1) || sm(1,1) != cplx(17, 0) ||
             sm(2,0) != cplx(14,2) || sm(2,1) != cplx(11, 1) ||
             sm(3,0) != cplx(15,3) || sm(3,1) != cplx(19, 2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) )\n"
                                        "( (18,1) (17, 0) )\n"
                                        "( (14,2) (11, 1) )\n"
                                        "( (15,3) (19, 2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(12,0) || herm(0,1) != cplx(18,-1) || herm(0,2) != cplx(14,-2) || herm(0,3) != cplx(15,-3) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(18,1) || herm(1,1) != cplx(17, 0) || herm(1,2) != cplx(11,-1) || herm(1,3) != cplx(19,-2) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,2) || herm(2,1) != cplx(11, 1) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(15,3) || herm(3,1) != cplx(19, 2) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )\n"
                                        "( (18,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )\n"
                                        "( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 2";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-15, 0);
         mat(0,2) = cplx(-13,-1);
         mat(0,3) = cplx(-15, 3);
         mat(1,0) = cplx(-13, 2);
         mat(1,1) = cplx(-13, 1);
         mat(1,2) = cplx( -6, 0);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18, 0) || sm(0,2) != cplx(14,2) || sm(0,3) != cplx(15,-3) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,-2) || sm(1,2) != cplx(11,0) || sm(1,3) != cplx(19, 1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18, 0) (14,2) (15,-3) )\n"
                                        "( (13,-2) (14,-2) (11,0) (19, 1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-13,-2);
         mat(1,0) = cplx(-15, 0);
         mat(1,1) = cplx(-13,-1);
         mat(2,0) = cplx(-13, 1);
         mat(2,1) = cplx( -6, 0);
         mat(3,0) = cplx(-15,-3);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13, 2) ||
             sm(1,0) != cplx(18, 0) || sm(1,1) != cplx(14, 2) ||
             sm(2,0) != cplx(14,-2) || sm(2,1) != cplx(11, 0) ||
             sm(3,0) != cplx(15, 3) || sm(3,1) != cplx(19,-1) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13, 2) )\n"
                                        "( (18, 0) (14, 2) )\n"
                                        "( (14,-2) (11, 0) )\n"
                                        "( (15, 3) (19,-1) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(12, 1) || herm(1,3) != cplx(13, 2) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx(12,-1) || herm(2,2) != cplx(18, 0) || herm(2,3) != cplx(14, 2) || herm(2,4) != cplx(15,-3) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx(13,-2) || herm(3,2) != cplx(14,-2) || herm(3,3) != cplx(11, 0) || herm(3,4) != cplx(19, 1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(15, 3) || herm(4,3) != cplx(19,-1) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )\n"
                                        "( (-2,-1) (13,-2) (14,-2) (11, 0) (19, 1) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 3";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-11,-2);
         mat(0,2) = cplx(-13, 0);
         mat(0,3) = cplx(-15,-1);
         mat(1,0) = cplx(-15, 2);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-15, 1);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(18,1) || sm(0,2) != cplx(14, 0) || sm(0,3) != cplx(11,1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14,0) || sm(1,2) != cplx(11,-1) || sm(1,3) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (18,1) (14, 0) (11,1) )\n"
                                        "( (13,-2) (14,0) (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-15,-2);
         mat(1,0) = cplx(-11, 2);
         mat(1,1) = cplx(-14, 0);
         mat(2,0) = cplx(-13, 0);
         mat(2,1) = cplx(-15,-1);
         mat(3,0) = cplx(-15, 1);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 30UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(13,2) ||
             sm(1,0) != cplx(18,-1) || sm(1,1) != cplx(14,0) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,1) ||
             sm(3,0) != cplx(11,-1) || sm(3,1) != cplx(19,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (13,2) )\n"
                                        "( (18,-1) (14,0) )\n"
                                        "( (14, 0) (11,1) )\n"
                                        "( (11,-1) (19,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2,1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0,0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1,1) || herm(2,4) != cplx(12, 1) || herm(2,5) != cplx(13, 2) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5,0) || herm(3,4) != cplx(18,-1) || herm(3,5) != cplx(14, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx(12,-1) || herm(4,3) != cplx(18,1) || herm(4,4) != cplx(14, 0) || herm(4,5) != cplx(11, 1) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(13,-2) || herm(5,3) != cplx(14,0) || herm(5,4) != cplx(11,-1) || herm(5,5) != cplx(19, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )\n"
                                        "( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )\n"
                                        "( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 4";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx( -5,-4);
         mat(0,1) = cplx(-18, 3);
         mat(0,2) = cplx(-11, 0);
         mat(0,3) = cplx(-10, 2);
         mat(1,0) = cplx(-15, 1);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-10,-2);
         mat(1,3) = cplx(-14, 0);
         mat(2,0) = cplx(-14,-3);
         mat(2,1) = cplx(-12,-4);
         mat(2,2) = cplx(-12,-1);
         mat(2,3) = cplx( -7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12, 1) || sm(0,1) != cplx(18,-3) || sm(0,2) != cplx(14,0) || sm(0,3) != cplx(11,-1) ||
             sm(1,0) != cplx(13,-2) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,1) || sm(1,3) != cplx(19, 0) ||
             sm(2,0) != cplx(19, 3) || sm(2,1) != cplx(11, 2) || sm(2,2) != cplx(12,1) || sm(2,3) != cplx(14,-4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12, 1) (18,-3) (14, 0) (11,-1) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx( -5, 4);
         mat(0,1) = cplx(-15,-1);
         mat(0,2) = cplx(-14, 3);
         mat(1,0) = cplx(-18,-3);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-12, 4);
         mat(2,0) = cplx(-11, 0);
         mat(2,1) = cplx(-10, 2);
         mat(2,2) = cplx(-12, 1);
         mat(3,0) = cplx(-10,-2);
         mat(3,1) = cplx(-14, 0);
         mat(3,2) = cplx( -7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );
         sm -= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 32UL );

         if( sm(0,0) != cplx(12,-1) || sm(0,1) != cplx(13, 2) || sm(0,2) != cplx(19,-3) ||
             sm(1,0) != cplx(18, 3) || sm(1,1) != cplx(14, 0) || sm(1,2) != cplx(11,-2) ||
             sm(2,0) != cplx(14, 0) || sm(2,1) != cplx(11,-1) || sm(2,2) != cplx(12,-1) ||
             sm(3,0) != cplx(11, 1) || sm(3,1) != cplx(19, 0) || sm(3,2) != cplx(14, 4) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (12,-1) (13, 2) (19,-3) )\n"
                                        "( (18, 3) (14, 0) (11,-2) )\n"
                                        "( (14, 0) (11,-1) (12,-1) )\n"
                                        "( (11, 1) (19, 0) (14, 4) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(12,-1) || herm(0,3) != cplx(13, 2) || herm(0,4) != cplx(19,-3) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx(18, 3) || herm(1,3) != cplx(14, 0) || herm(1,4) != cplx(11,-2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(12, 1) || herm(2,1) != cplx(18,-3) || herm(2,2) != cplx(14, 0) || herm(2,3) != cplx(11,-1) || herm(2,4) != cplx(12,-1) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(13,-2) || herm(3,1) != cplx(14, 0) || herm(3,2) != cplx(11, 1) || herm(3,3) != cplx(19, 0) || herm(3,4) != cplx(14, 4) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(19, 3) || herm(4,1) != cplx(11, 2) || herm(4,2) != cplx(12, 1) || herm(4,3) != cplx(14,-4) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )\n"
                                        "( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )\n"
                                        "( (13,-2) (14, 0) (11, 1) (19, 0) (14, 4) ( 0, 0) )\n"
                                        "( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (22,1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 5";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-11,0);
         mat(0,1) = cplx(-22,0);
         mat(0,2) = cplx( -7,5);
         mat(0,3) = cplx(-17,4);
         mat(1,0) = cplx(-26,0);
         mat(1,1) = cplx(-15,0);
         mat(1,2) = cplx(-11,1);
         mat(1,3) = cplx(-19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-11, 0);
         mat(0,1) = cplx(-22, 0);
         mat(1,0) = cplx(-26, 0);
         mat(1,1) = cplx(-15, 0);
         mat(2,0) = cplx( -7,-5);
         mat(2,1) = cplx(-11,-1);
         mat(3,0) = cplx(-17,-4);
         mat(3,1) = cplx(-19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (22,-2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 6";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-15, 0);
         mat(0,2) = cplx(-13,-1);
         mat(0,3) = cplx(-15, 3);
         mat(1,0) = cplx(-13, 2);
         mat(1,1) = cplx(-21, 1);
         mat(1,2) = cplx( -6, 0);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-13,-2);
         mat(1,0) = cplx(-15, 0);
         mat(1,1) = cplx(-13,-1);
         mat(2,0) = cplx(-21, 1);
         mat(2,1) = cplx( -6, 0);
         mat(3,0) = cplx(-15,-3);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (22,-1) (19, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 7";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-11,-2);
         mat(0,2) = cplx(-13, 0);
         mat(0,3) = cplx(-15,-1);
         mat(1,0) = cplx(-15, 2);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-26, 1);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-15,-2);
         mat(1,0) = cplx(-11, 2);
         mat(1,1) = cplx(-14, 0);
         mat(2,0) = cplx(-13, 0);
         mat(2,1) = cplx(-15,-1);
         mat(3,0) = cplx(-26, 1);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (22, 1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 8";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx( -5,-4);
         mat(0,1) = cplx(-18, 3);
         mat(0,2) = cplx(-11, 0);
         mat(0,3) = cplx(-10, 2);
         mat(1,0) = cplx(-15, 1);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-21,-2);
         mat(1,3) = cplx(-14, 0);
         mat(2,0) = cplx(-14,-3);
         mat(2,1) = cplx(-12,-4);
         mat(2,2) = cplx(-12,-1);
         mat(2,3) = cplx( -7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx( -5, 4);
         mat(0,1) = cplx(-15,-1);
         mat(0,2) = cplx(-14, 3);
         mat(1,0) = cplx(-18,-3);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-12, 4);
         mat(2,0) = cplx(-11, 0);
         mat(2,1) = cplx(-10, 2);
         mat(2,2) = cplx(-12, 1);
         mat(3,0) = cplx(-21,-2);
         mat(3,1) = cplx(-14, 0);
         mat(3,2) = cplx( -7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12, 0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,-1) (17, 0) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14, 2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15, 3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 9";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-11,0);
         mat(0,1) = cplx(-22,0);
         mat(0,2) = cplx( -7,5);
         mat(0,3) = cplx(-17,4);
         mat(1,0) = cplx(-22,2);
         mat(1,1) = cplx(-15,0);
         mat(1,2) = cplx(-11,1);
         mat(1,3) = cplx(-19,2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-11, 0);
         mat(0,1) = cplx(-22, 0);
         mat(1,0) = cplx(-22, 2);
         mat(1,1) = cplx(-15, 0);
         mat(2,0) = cplx( -7,-5);
         mat(2,1) = cplx(-11,-1);
         mat(3,0) = cplx(-17,-4);
         mat(3,1) = cplx(-19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12,1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18,0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,2) (11, 0) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15,3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 10";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-15, 0);
         mat(0,2) = cplx(-13,-1);
         mat(0,3) = cplx(-15, 3);
         mat(1,0) = cplx(-13, 2);
         mat(1,1) = cplx(-13, 3);
         mat(1,2) = cplx( -6, 0);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-13,-2);
         mat(1,0) = cplx(-15, 0);
         mat(1,1) = cplx(-13,-1);
         mat(2,0) = cplx(-13, 3);
         mat(2,1) = cplx( -6, 0);
         mat(3,0) = cplx(-15,-3);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11, 1) (19, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 11";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-11,-2);
         mat(0,2) = cplx(-13, 0);
         mat(0,3) = cplx(-15,-1);
         mat(1,0) = cplx(-15, 2);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-15,-1);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-15,-2);
         mat(1,0) = cplx(-11, 2);
         mat(1,1) = cplx(-14, 0);
         mat(2,0) = cplx(-13, 0);
         mat(2,1) = cplx(-15,-1);
         mat(3,0) = cplx(-15,-1);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11,-1) (19, 0) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 12";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx( -5,-4);
         mat(0,1) = cplx(-18, 3);
         mat(0,2) = cplx(-11, 0);
         mat(0,3) = cplx(-10, 2);
         mat(1,0) = cplx(-15, 1);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-10, 0);
         mat(1,3) = cplx(-14, 0);
         mat(2,0) = cplx(-14,-3);
         mat(2,1) = cplx(-12,-4);
         mat(2,2) = cplx(-12,-1);
         mat(2,3) = cplx( -7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx( -5, 4);
         mat(0,1) = cplx(-15,-1);
         mat(0,2) = cplx(-14, 3);
         mat(1,0) = cplx(-18,-3);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-12, 4);
         mat(2,0) = cplx(-11, 0);
         mat(2,1) = cplx(-10, 2);
         mat(2,2) = cplx(-12, 1);
         mat(3,0) = cplx(-10, 0);
         mat(3,1) = cplx(-14, 0);
         mat(3,2) = cplx( -7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (12,0) (18,-1) (14,-2) (15,-3) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (18,1) (17, 1) (11,-1) (19,-2) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,2) (11, 1) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (15,3) (19, 2) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 13";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-11, 0);
         mat(0,1) = cplx(-22, 0);
         mat(0,2) = cplx( -7, 5);
         mat(0,3) = cplx(-17, 4);
         mat(1,0) = cplx(-22, 0);
         mat(1,1) = cplx(-15,-1);
         mat(1,2) = cplx(-11, 1);
         mat(1,3) = cplx(-19, 2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-11, 0);
         mat(0,1) = cplx(-22, 0);
         mat(1,0) = cplx(-22, 0);
         mat(1,1) = cplx(-15,-1);
         mat(2,0) = cplx( -7,-5);
         mat(2,1) = cplx(-11,-1);
         mat(3,0) = cplx(-17,-4);
         mat(3,1) = cplx(-19,-2);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (12, 1) (13, 2) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) (12,-1) (18, 0) (14, 2) (15,-3) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) (13,-2) (14,-2) (11, 1) (19, 1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (15, 3) (19,-1) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 14";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-15, 0);
         mat(0,2) = cplx(-13,-1);
         mat(0,3) = cplx(-15, 3);
         mat(1,0) = cplx(-13, 2);
         mat(1,1) = cplx(-13, 1);
         mat(1,2) = cplx( -6,-1);
         mat(1,3) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-13,-2);
         mat(1,0) = cplx(-15, 0);
         mat(1,1) = cplx(-13,-1);
         mat(2,0) = cplx(-13, 1);
         mat(2,1) = cplx( -6,-1);
         mat(3,0) = cplx(-15,-3);
         mat(3,1) = cplx(-12, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2,1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0,0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1,1) (12, 1) (13, 2) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5,0) (18,-1) (14, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) (12,-1) (18,1) (14, 0) (11, 1) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (13,-2) (14,0) (11,-1) (19, 1) )
   {
      test_ = "Sparse matrix subtraction assignment test 15";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(-12, 1);
         mat(0,1) = cplx(-11,-2);
         mat(0,2) = cplx(-13, 0);
         mat(0,3) = cplx(-15,-1);
         mat(1,0) = cplx(-15, 2);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-15, 1);
         mat(1,3) = cplx(-12,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(-12,-1);
         mat(0,1) = cplx(-15,-2);
         mat(1,0) = cplx(-11, 2);
         mat(1,1) = cplx(-14, 0);
         mat(2,0) = cplx(-13, 0);
         mat(2,1) = cplx(-15,-1);
         mat(3,0) = cplx(-15, 1);
         mat(3,1) = cplx(-12,-1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (12,-1) (13, 2) (19,-3) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) (18, 3) (14, 0) (11,-2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (12, 1) (18,-3) (14, 0) (11,-1) (12,-1) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (13,-2) (14, 0) (11, 1) (19, 1) (14, 4) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (19, 3) (11, 2) (12, 1) (14,-4) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix subtraction assignment test 16";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx( -5,-4);
         mat(0,1) = cplx(-18, 3);
         mat(0,2) = cplx(-11, 0);
         mat(0,3) = cplx(-10, 2);
         mat(1,0) = cplx(-15, 1);
         mat(1,1) = cplx(-14, 0);
         mat(1,2) = cplx(-10,-2);
         mat(1,3) = cplx(-14,-1);
         mat(2,0) = cplx(-14,-3);
         mat(2,1) = cplx(-12,-4);
         mat(2,2) = cplx(-12,-1);
         mat(2,3) = cplx( -7, 3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx( -5, 4);
         mat(0,1) = cplx(-15,-1);
         mat(0,2) = cplx(-14, 3);
         mat(1,0) = cplx(-18,-3);
         mat(1,1) = cplx(-14,-1);
         mat(1,2) = cplx(-12, 4);
         mat(2,0) = cplx(-11, 0);
         mat(2,1) = cplx(-10, 2);
         mat(2,2) = cplx(-12, 1);
         mat(3,0) = cplx(-10,-2);
         mat(3,1) = cplx(-14,-1);
         mat(3,2) = cplx( -7,-3);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm -= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the Schur product assignment to a submatrix of a HermitianMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the Schur product assignment to a submatrix of a HermitianMatrix.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename HT >  // Type of the Hermitian matrix
void SubmatrixComplexTest::testSchurAssign()
{
   //=====================================================================================
   // Dense matrix Schur product assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (11,  0) (20, 5) (28,12) (16,-8) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (20, -5) (12, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (28,-12) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (16,  8) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,  0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,  0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix Schur product assignment test 1";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-5, 0);
         mat(0,2) = cplx( 4, 0);
         mat(0,3) = cplx(-8, 0);
         mat(1,0) = cplx(-5, 0);
         mat(1,1) = cplx( 6, 0);
         mat(1,2) = cplx(99,99);
         mat(1,3) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx(11, 0) || sm(0,1) != cplx(20,5) || sm(0,2) != cplx(28,12) || sm(0,3) != cplx(16,-8) ||
             sm(1,0) != cplx(20,-5) || sm(1,1) != cplx(12,0) || sm(1,2) != cplx( 0, 0) || sm(1,3) != cplx( 0, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (11, 0) (20,5) (28,12) (16,-8) )\n"
                                        "( (20,-5) (12,0) ( 0, 0) ( 0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(11,  0) || herm(0,1) != cplx(20, 5) || herm(0,2) != cplx(28,12) || herm(0,3) != cplx(16,-8) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(20, -5) || herm(1,1) != cplx(12, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(28,-12) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(16,  8) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,  0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,  0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (11,  0) (20, 5) (28,12) (16,-8) ( 5,0) ( 0, 0) )\n"
                                        "( (20, -5) (12, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )\n"
                                        "( (28,-12) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (16,  8) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,  0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,  0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-5, 0);
         mat(1,0) = cplx(-5, 0);
         mat(1,1) = cplx( 6, 0);
         mat(2,0) = cplx( 4, 0);
         mat(2,1) = cplx(99,99);
         mat(3,0) = cplx(-8, 0);
         mat(3,1) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx(11,  0) || sm(0,1) != cplx(20,5) ||
             sm(1,0) != cplx(20, -5) || sm(1,1) != cplx(12,0) ||
             sm(2,0) != cplx(28,-12) || sm(2,1) != cplx( 0,0) ||
             sm(3,0) != cplx(16,  8) || sm(3,1) != cplx( 0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (11,  0) (20,5) )\n"
                                        "( (20, -5) (12,0) )\n"
                                        "( (28,-12) ( 0,0) )\n"
                                        "( (16,  8) ( 0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(11,  0) || herm(0,1) != cplx(20, 5) || herm(0,2) != cplx(28,12) || herm(0,3) != cplx(16,-8) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(20, -5) || herm(1,1) != cplx(12, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(28,-12) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(16,  8) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,  0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,  0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (11,  0) (20, 5) (28,12) (16,-8) ( 5,0) ( 0, 0) )\n"
                                        "( (20, -5) (12, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )\n"
                                        "( (28,-12) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (16,  8) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,  0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,  0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,  3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) (18,  0) (14,14) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) (14,-14) (20, 0) (21, 3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0,  0) (21,-3) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix Schur product assignment test 2";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 6, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(99,99);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx( 4, 0);
         mat(1,3) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx(0,0) || sm(0,1) != cplx(18,  0) || sm(0,2) != cplx(14,14) || sm(0,3) != cplx( 0,0) ||
             sm(1,0) != cplx(0,0) || sm(1,1) != cplx(14,-14) || sm(1,2) != cplx(20, 0) || sm(1,3) != cplx(21,3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (0,0) (18,  0) (14,14) ( 0,0) )\n"
                                        "( (0,0) (14,-14) (20, 0) (21,3) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7,  3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0,  0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(18,  0) || herm(2,3) != cplx(14,14) || herm(2,4) != cplx( 0, 0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx(14,-14) || herm(3,3) != cplx(20, 0) || herm(3,4) != cplx(21, 3) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0,  0) || herm(4,3) != cplx(21,-3) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2,  0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7,  3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) (18,  0) (14,14) ( 0, 0) (-2, 0) )\n"
                                        "( (-2,-1) ( 0, 0) (14,-14) (20, 0) (21, 3) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) ( 0,  0) (21,-3) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(99,99);
         mat(1,0) = cplx( 6, 0);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx( 4, 0);
         mat(3,0) = cplx(99,99);
         mat(3,1) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx( 0,  0) || sm(0,1) != cplx( 0, 0) ||
             sm(1,0) != cplx(18,  0) || sm(1,1) != cplx(14,14) ||
             sm(2,0) != cplx(14,-14) || sm(2,1) != cplx(20, 0) ||
             sm(3,0) != cplx( 0,  0) || sm(3,1) != cplx(21,-3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 0,  0) ( 0, 0) )\n"
                                        "( (18,  0) (14,14) )\n"
                                        "( (14,-14) (20, 0) )\n"
                                        "( ( 0,  0) (21,-3) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7,  3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0,  0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(18,  0) || herm(2,3) != cplx(14,14) || herm(2,4) != cplx( 0, 0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx(14,-14) || herm(3,3) != cplx(20, 0) || herm(3,4) != cplx(21, 3) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0,  0) || herm(4,3) != cplx(21,-3) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2,  0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7,  3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) (18,  0) (14,14) ( 0, 0) (-2, 0) )\n"
                                        "( (-2,-1) ( 0, 0) (14,-14) (20, 0) (21, 3) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) ( 0,  0) (21,-3) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (16, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) (21,3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) (21,-3) (14,0) (20, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (16, 0) ( 0, 0) (20,0) (28, 0) )
   {
      test_ = "Dense matrix Schur product assignment test 3";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 3, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(-5, 0);
         mat(1,0) = cplx(-8, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-5, 0);
         mat(1,3) = cplx( 4, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx( 0,0) || sm(0,1) != cplx(21,-3) || sm(0,2) != cplx(14,0) || sm(0,3) != cplx(20,0) ||
             sm(1,0) != cplx(16,0) || sm(1,1) != cplx( 0, 0) || sm(1,2) != cplx(20,0) || sm(1,3) != cplx(28,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 0,0) (21,-3) (14,0) (20,0) )\n"
                                        "( (16,0) ( 0, 0) (20,0) (28,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(16, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx(21,3) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx(21,-3) || herm(4,4) != cplx(14,0) || herm(4,5) != cplx(20, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(16, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(20,0) || herm(5,5) != cplx(28, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0, 0) (16, 0) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) (21, 3) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) ( 0, 0) (21,-3) (14, 0) (20, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (16, 0) ( 0, 0) (20, 0) (28, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(-8, 0);
         mat(1,0) = cplx( 3, 0);
         mat(1,1) = cplx(99,99);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(-5, 0);
         mat(3,0) = cplx(-5, 0);
         mat(3,1) = cplx( 4, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx( 0,0) || sm(0,1) != cplx(16,0) ||
             sm(1,0) != cplx(21,3) || sm(1,1) != cplx( 0,0) ||
             sm(2,0) != cplx(14,0) || sm(2,1) != cplx(20,0) ||
             sm(3,0) != cplx(20,0) || sm(3,1) != cplx(28,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 0,0) (16,0) )\n"
                                        "( (21,3) ( 0,0) )\n"
                                        "( (14,0) (20,0) )\n"
                                        "( (20,0) (28,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(16, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx(21,3) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx(21,-3) || herm(4,4) != cplx(14,0) || herm(4,5) != cplx(20, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(16, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(20,0) || herm(5,5) != cplx(28, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0, 0) (16, 0) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) (21, 3) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) ( 0, 0) (21,-3) (14, 0) (20, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (16, 0) ( 0, 0) (20, 0) (28, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (14,  6) (18,-9) (25, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) ( 7,14) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,-6) ( 0, 0) (18,  0) (11,11) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (18, 9) ( 0, 0) (11,-11) (20, 0) (14, 2) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (25, 0) ( 7,14) ( 0,  0) (14,-2) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix Schur product assignment test 4";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(99,99);
         mat(0,2) = cplx( 6, 0);
         mat(0,3) = cplx(11, 0);
         mat(1,0) = cplx(-9, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(11, 0);
         mat(1,3) = cplx( 4, 0);
         mat(2,0) = cplx( 5, 0);
         mat(2,1) = cplx(-7, 0);
         mat(2,2) = cplx(99,99);
         mat(2,3) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx(14,-6) || sm(0,1) != cplx(0, 0) || sm(0,2) != cplx(18,  0) || sm(0,3) != cplx(11,11) ||
             sm(1,0) != cplx(18, 9) || sm(1,1) != cplx(0, 0) || sm(1,2) != cplx(11,-11) || sm(1,3) != cplx(20, 0) ||
             sm(2,0) != cplx(25, 0) || sm(2,1) != cplx(7,14) || sm(2,2) != cplx( 0,  0) || sm(2,3) != cplx(14,-2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (14,-6) (0, 0) (18,  0) (11,11) )\n"
                                        "( (18, 9) (0, 0) (11,-11) (20, 0) )\n"
                                        "( (25, 0) (7,14) ( 0,  0) (14,-2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(14,  6) || herm(0,3) != cplx(18,-9) || herm(0,4) != cplx(25,  0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0,  0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx( 7,-14) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(18,  0) || herm(2,3) != cplx(11,11) || herm(2,4) != cplx( 0,  0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(18, 9) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx(11,-11) || herm(3,3) != cplx(20, 0) || herm(3,4) != cplx(14,  2) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(25, 0) || herm(4,1) != cplx( 7,14) || herm(4,2) != cplx( 0,  0) || herm(4,3) != cplx(14,-2) || herm(4,4) != cplx( 1,  0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2,  0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,  0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (14,  6) (18,-9) (25,  0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) ( 7,-14) ( 8,-2) )\n"
                                        "( (14,-6) ( 0, 0) (18,  0) (11,11) ( 0,  0) (-2, 0) )\n"
                                        "( (18, 9) ( 0, 0) (11,-11) (20, 0) (14,  2) ( 0, 0) )\n"
                                        "( (25, 0) ( 7,14) ( 0,  0) (14,-2) ( 1,  0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4,  0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(-9, 0);
         mat(0,2) = cplx( 5, 0);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-7, 0);
         mat(2,0) = cplx( 6, 0);
         mat(2,1) = cplx(11, 0);
         mat(2,2) = cplx(99,99);
         mat(3,0) = cplx(11, 0);
         mat(3,1) = cplx( 4, 0);
         mat(3,2) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx(14,  6) || sm(0,1) != cplx(18,-9) || sm(0,2) != cplx(25,  0) ||
             sm(1,0) != cplx( 0,  0) || sm(1,1) != cplx( 0, 0) || sm(1,2) != cplx( 7,-14) ||
             sm(2,0) != cplx(18,  0) || sm(2,1) != cplx(11,11) || sm(2,2) != cplx( 0,  0) ||
             sm(3,0) != cplx(11,-11) || sm(3,1) != cplx(20, 0) || sm(3,2) != cplx(14,  2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (14,  6) (18,-9) (25,  0) )\n"
                                        "( ( 0,  0) ( 0, 0) ( 7,-14) )\n"
                                        "( (18,  0) (11,11) ( 0,  0) )\n"
                                        "( (11,-11) (20, 0) (14,  2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(14,  6) || herm(0,3) != cplx(18,-9) || herm(0,4) != cplx(25,  0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0,  0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx( 7,-14) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(18,  0) || herm(2,3) != cplx(11,11) || herm(2,4) != cplx( 0,  0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(18, 9) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx(11,-11) || herm(3,3) != cplx(20, 0) || herm(3,4) != cplx(14,  2) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(25, 0) || herm(4,1) != cplx( 7,14) || herm(4,2) != cplx( 0,  0) || herm(4,3) != cplx(14,-2) || herm(4,4) != cplx( 1,  0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2,  0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,  0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (14,  6) (18,-9) (25,  0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) ( 7,-14) ( 8,-2) )\n"
                                        "( (14,-6) ( 0, 0) (18,  0) (11,11) ( 0,  0) (-2, 0) )\n"
                                        "( (18, 9) ( 0, 0) (11,-11) (20, 0) (14,  2) ( 0, 0) )\n"
                                        "( (25, 0) ( 7,14) ( 0,  0) (14,-2) ( 1,  0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4,  0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (11,  0) (20, 5) (28,12) (16,-8) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (24, -6) (12, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (28,-12) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (16,  8) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,  0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,  0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix Schur product assignment test 5";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-5, 0);
         mat(0,2) = cplx( 4, 0);
         mat(0,3) = cplx(-8, 0);
         mat(1,0) = cplx(-6, 0);
         mat(1,1) = cplx( 6, 0);
         mat(1,2) = cplx(99,99);
         mat(1,3) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-6, 0);
         mat(1,0) = cplx(-5, 0);
         mat(1,1) = cplx( 6, 0);
         mat(2,0) = cplx( 4, 0);
         mat(2,1) = cplx(99,99);
         mat(3,0) = cplx(-8, 0);
         mat(3,1) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,  3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) (18,  0) (14,14) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) (22,-22) (20, 0) (21, 3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0,  0) (21,-3) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix Schur product assignment test 6";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 6, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(99,99);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(22, 0);
         mat(1,2) = cplx( 4, 0);
         mat(1,3) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(99,99);
         mat(1,0) = cplx( 6, 0);
         mat(1,1) = cplx(22, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx( 4, 0);
         mat(3,0) = cplx(99,99);
         mat(3,1) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (16, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) (21,3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) (21,-3) (14,0) (20, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (16, 0) ( 0, 0) (24,0) (28, 0) )
   {
      test_ = "Dense matrix Schur product assignment test 7";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 3, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(-5, 0);
         mat(1,0) = cplx(-8, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-6, 0);
         mat(1,3) = cplx( 4, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(-8, 0);
         mat(1,0) = cplx( 3, 0);
         mat(1,1) = cplx(99,99);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(-6, 0);
         mat(3,0) = cplx(-5, 0);
         mat(3,1) = cplx( 4, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (14,  6) (18,-9) (25, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) ( 7,14) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,-6) ( 0, 0) (18,  0) (11,11) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (18, 9) ( 0, 0) (22,-22) (20, 0) (14, 2) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (25, 0) ( 7,14) ( 0,  0) (14,-2) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix Schur product assignment test 8";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(99,99);
         mat(0,2) = cplx( 6, 0);
         mat(0,3) = cplx(11, 0);
         mat(1,0) = cplx(-9, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(22, 0);
         mat(1,3) = cplx( 4, 0);
         mat(2,0) = cplx( 5, 0);
         mat(2,1) = cplx(-7, 0);
         mat(2,2) = cplx(99,99);
         mat(2,3) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(-9, 0);
         mat(0,2) = cplx( 5, 0);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-7, 0);
         mat(2,0) = cplx( 6, 0);
         mat(2,1) = cplx(22, 0);
         mat(2,2) = cplx(99,99);
         mat(3,0) = cplx(11, 0);
         mat(3,1) = cplx( 4, 0);
         mat(3,2) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (11,  0) (20, 5) (28,12) (16,-8) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (20, -5) (12, 4) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (28,-12) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (16,  8) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,  0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,  0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Dense matrix Schur product assignment test 9";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-5, 0);
         mat(0,2) = cplx( 4, 0);
         mat(0,3) = cplx(-8, 0);
         mat(1,0) = cplx(-5, 0);
         mat(1,1) = cplx( 6, 2);
         mat(1,2) = cplx(99,99);
         mat(1,3) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-5, 0);
         mat(1,0) = cplx(-5, 0);
         mat(1,1) = cplx( 6, 2);
         mat(2,0) = cplx( 4, 0);
         mat(2,1) = cplx(99,99);
         mat(3,0) = cplx(-8, 0);
         mat(3,1) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,  3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) (18,  0) (14,14) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) (14,-14) (20, 5) (21, 3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0,  0) (21,-3) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix Schur product assignment test 10";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 6, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(99,99);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx( 4, 1);
         mat(1,3) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(99,99);
         mat(1,0) = cplx( 6, 0);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx( 4, 1);
         mat(3,0) = cplx(99,99);
         mat(3,1) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (16, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) (21,3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) (21,-3) (14,0) (20, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (16, 0) ( 0, 0) (20,0) (28, 7) )
   {
      test_ = "Dense matrix Schur product assignment test 11";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 3, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(-5, 0);
         mat(1,0) = cplx(-8, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-5, 0);
         mat(1,3) = cplx( 4, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(-8, 0);
         mat(1,0) = cplx( 3, 0);
         mat(1,1) = cplx(99,99);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(-5, 0);
         mat(3,0) = cplx(-5, 0);
         mat(3,1) = cplx( 4, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (14,  6) (18,-9) (25, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) ( 7,14) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,-6) ( 0, 0) (18,  0) (11,11) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (18, 9) ( 0, 0) (11,-11) (20, 5) (14, 2) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (25, 0) ( 7,14) ( 0,  0) (14,-2) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Dense matrix Schur product assignment test 12";

      {
         blaze::DynamicMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(99,99);
         mat(0,2) = cplx( 6, 0);
         mat(0,3) = cplx(11, 0);
         mat(1,0) = cplx(-9, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(11, 0);
         mat(1,3) = cplx( 4, 1);
         mat(2,0) = cplx( 5, 0);
         mat(2,1) = cplx(-7, 0);
         mat(2,2) = cplx(99,99);
         mat(2,3) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::DynamicMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(-9, 0);
         mat(0,2) = cplx( 5, 0);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-7, 0);
         mat(2,0) = cplx( 6, 0);
         mat(2,1) = cplx(11, 0);
         mat(2,2) = cplx(99,99);
         mat(3,0) = cplx(11, 0);
         mat(3,1) = cplx( 4, 1);
         mat(3,2) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }


   //=====================================================================================
   // Sparse matrix Schur product assignment
   //=====================================================================================

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (11,  0) (20, 5) (28,12) (16,-8) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (20, -5) (12, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (28,-12) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (16,  8) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,  0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,  0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix Schur product assignment test 1";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-5, 0);
         mat(0,2) = cplx( 4, 0);
         mat(0,3) = cplx(-8, 0);
         mat(1,0) = cplx(-5, 0);
         mat(1,1) = cplx( 6, 0);
         mat(1,2) = cplx(99,99);
         mat(1,3) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx(11, 0) || sm(0,1) != cplx(20,5) || sm(0,2) != cplx(28,12) || sm(0,3) != cplx(16,-8) ||
             sm(1,0) != cplx(20,-5) || sm(1,1) != cplx(12,0) || sm(1,2) != cplx( 0, 0) || sm(1,3) != cplx( 0, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (11, 0) (20,5) (28,12) (16,-8) )\n"
                                        "( (20,-5) (12,0) ( 0, 0) ( 0, 0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(11,  0) || herm(0,1) != cplx(20, 5) || herm(0,2) != cplx(28,12) || herm(0,3) != cplx(16,-8) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(20, -5) || herm(1,1) != cplx(12, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(28,-12) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(16,  8) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,  0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,  0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (11,  0) (20, 5) (28,12) (16,-8) ( 5,0) ( 0, 0) )\n"
                                        "( (20, -5) (12, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )\n"
                                        "( (28,-12) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (16,  8) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,  0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,  0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-5, 0);
         mat(1,0) = cplx(-5, 0);
         mat(1,1) = cplx( 6, 0);
         mat(2,0) = cplx( 4, 0);
         mat(2,1) = cplx(99,99);
         mat(3,0) = cplx(-8, 0);
         mat(3,1) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx(11,  0) || sm(0,1) != cplx(20,5) ||
             sm(1,0) != cplx(20, -5) || sm(1,1) != cplx(12,0) ||
             sm(2,0) != cplx(28,-12) || sm(2,1) != cplx( 0,0) ||
             sm(3,0) != cplx(16,  8) || sm(3,1) != cplx( 0,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (11,  0) (20,5) )\n"
                                        "( (20, -5) (12,0) )\n"
                                        "( (28,-12) ( 0,0) )\n"
                                        "( (16,  8) ( 0,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx(11,  0) || herm(0,1) != cplx(20, 5) || herm(0,2) != cplx(28,12) || herm(0,3) != cplx(16,-8) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(20, -5) || herm(1,1) != cplx(12, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(28,-12) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(16,  8) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx( 7,1) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5,  0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx( 7,-1) || herm(4,4) != cplx( 1,0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0,  0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( (11,  0) (20, 5) (28,12) (16,-8) ( 5,0) ( 0, 0) )\n"
                                        "( (20, -5) (12, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )\n"
                                        "( (28,-12) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )\n"
                                        "( (16,  8) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )\n"
                                        "( ( 5,  0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )\n"
                                        "( ( 0,  0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,  3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) (18,  0) (14,14) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) (14,-14) (20, 0) (21, 3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0,  0) (21,-3) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix Schur product assignment test 2";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 6, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(99,99);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx( 4, 0);
         mat(1,3) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx(0,0) || sm(0,1) != cplx(18,  0) || sm(0,2) != cplx(14,14) || sm(0,3) != cplx( 0,0) ||
             sm(1,0) != cplx(0,0) || sm(1,1) != cplx(14,-14) || sm(1,2) != cplx(20, 0) || sm(1,3) != cplx(21,3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (0,0) (18,  0) (14,14) ( 0,0) )\n"
                                        "( (0,0) (14,-14) (20, 0) (21,3) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7,  3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0,  0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(18,  0) || herm(2,3) != cplx(14,14) || herm(2,4) != cplx( 0, 0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx(14,-14) || herm(3,3) != cplx(20, 0) || herm(3,4) != cplx(21, 3) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0,  0) || herm(4,3) != cplx(21,-3) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2,  0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7,  3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) (18,  0) (14,14) ( 0, 0) (-2, 0) )\n"
                                        "( (-2,-1) ( 0, 0) (14,-14) (20, 0) (21, 3) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) ( 0,  0) (21,-3) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(99,99);
         mat(1,0) = cplx( 6, 0);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx( 4, 0);
         mat(3,0) = cplx(99,99);
         mat(3,1) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx( 0,  0) || sm(0,1) != cplx( 0, 0) ||
             sm(1,0) != cplx(18,  0) || sm(1,1) != cplx(14,14) ||
             sm(2,0) != cplx(14,-14) || sm(2,1) != cplx(20, 0) ||
             sm(3,0) != cplx( 0,  0) || sm(3,1) != cplx(21,-3) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 0,  0) ( 0, 0) )\n"
                                        "( (18,  0) (14,14) )\n"
                                        "( (14,-14) (20, 0) )\n"
                                        "( ( 0,  0) (21,-3) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7,  3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5, 0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0,  0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1, 2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(18,  0) || herm(2,3) != cplx(14,14) || herm(2,4) != cplx( 0, 0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx(14,-14) || herm(3,3) != cplx(20, 0) || herm(3,4) != cplx(21, 3) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0,  0) || herm(4,3) != cplx(21,-3) || herm(4,4) != cplx( 1, 0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2,  0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4, 0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7,  3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) (18,  0) (14,14) ( 0, 0) (-2, 0) )\n"
                                        "( (-2,-1) ( 0, 0) (14,-14) (20, 0) (21, 3) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) ( 0,  0) (21,-3) ( 1, 0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (16, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) (21,3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) (21,-3) (14,0) (20, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (16, 0) ( 0, 0) (20,0) (28, 0) )
   {
      test_ = "Sparse matrix Schur product assignment test 3";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 3, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(-5, 0);
         mat(1,0) = cplx(-8, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-5, 0);
         mat(1,3) = cplx( 4, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx( 0,0) || sm(0,1) != cplx(21,-3) || sm(0,2) != cplx(14,0) || sm(0,3) != cplx(20,0) ||
             sm(1,0) != cplx(16,0) || sm(1,1) != cplx( 0, 0) || sm(1,2) != cplx(20,0) || sm(1,3) != cplx(28,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 0,0) (21,-3) (14,0) (20,0) )\n"
                                        "( (16,0) ( 0, 0) (20,0) (28,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(16, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx(21,3) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx(21,-3) || herm(4,4) != cplx(14,0) || herm(4,5) != cplx(20, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(16, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(20,0) || herm(5,5) != cplx(28, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0, 0) (16, 0) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) (21, 3) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) ( 0, 0) (21,-3) (14, 0) (20, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (16, 0) ( 0, 0) (20, 0) (28, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(-8, 0);
         mat(1,0) = cplx( 3, 0);
         mat(1,1) = cplx(99,99);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(-5, 0);
         mat(3,0) = cplx(-5, 0);
         mat(3,1) = cplx( 4, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx( 0,0) || sm(0,1) != cplx(16,0) ||
             sm(1,0) != cplx(21,3) || sm(1,1) != cplx( 0,0) ||
             sm(2,0) != cplx(14,0) || sm(2,1) != cplx(20,0) ||
             sm(3,0) != cplx(20,0) || sm(3,1) != cplx(28,0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( ( 0,0) (16,0) )\n"
                                        "( (21,3) ( 0,0) )\n"
                                        "( (14,0) (20,0) )\n"
                                        "( (20,0) (28,0) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx( 7, 3) || herm(0,3) != cplx(-2, 1) || herm(0,4) != cplx( 5,0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0, 0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx(-1,2) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx( 7,-3) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx( 3, 0) || herm(2,3) != cplx( 1, 1) || herm(2,4) != cplx( 0,0) || herm(2,5) != cplx(16, 0) ||
             herm(3,0) != cplx(-2,-1) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx( 1,-1) || herm(3,3) != cplx( 5, 0) || herm(3,4) != cplx(21,3) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx( 5, 0) || herm(4,1) != cplx(-1,-2) || herm(4,2) != cplx( 0, 0) || herm(4,3) != cplx(21,-3) || herm(4,4) != cplx(14,0) || herm(4,5) != cplx(20, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(16, 0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(20,0) || herm(5,5) != cplx(28, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5, 0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1, 2) ( 8,-2) )\n"
                                        "( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0, 0) (16, 0) )\n"
                                        "( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) (21, 3) ( 0, 0) )\n"
                                        "( ( 5, 0) (-1,-2) ( 0, 0) (21,-3) (14, 0) (20, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (16, 0) ( 0, 0) (20, 0) (28, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (14,  6) (18,-9) (25, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) ( 7,14) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,-6) ( 0, 0) (18,  0) (11,11) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (18, 9) ( 0, 0) (11,-11) (20, 0) (14, 2) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (25, 0) ( 7,14) ( 0,  0) (14,-2) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix Schur product assignment test 4";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(99,99);
         mat(0,2) = cplx( 6, 0);
         mat(0,3) = cplx(11, 0);
         mat(1,0) = cplx(-9, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(11, 0);
         mat(1,3) = cplx( 4, 0);
         mat(2,0) = cplx( 5, 0);
         mat(2,1) = cplx(-7, 0);
         mat(2,2) = cplx(99,99);
         mat(2,3) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx(14,-6) || sm(0,1) != cplx(0, 0) || sm(0,2) != cplx(18,  0) || sm(0,3) != cplx(11,11) ||
             sm(1,0) != cplx(18, 9) || sm(1,1) != cplx(0, 0) || sm(1,2) != cplx(11,-11) || sm(1,3) != cplx(20, 0) ||
             sm(2,0) != cplx(25, 0) || sm(2,1) != cplx(7,14) || sm(2,2) != cplx( 0,  0) || sm(2,3) != cplx(14,-2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (14,-6) (0, 0) (18,  0) (11,11) )\n"
                                        "( (18, 9) (0, 0) (11,-11) (20, 0) )\n"
                                        "( (25, 0) (7,14) ( 0,  0) (14,-2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(14,  6) || herm(0,3) != cplx(18,-9) || herm(0,4) != cplx(25,  0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0,  0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx( 7,-14) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(18,  0) || herm(2,3) != cplx(11,11) || herm(2,4) != cplx( 0,  0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(18, 9) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx(11,-11) || herm(3,3) != cplx(20, 0) || herm(3,4) != cplx(14,  2) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(25, 0) || herm(4,1) != cplx( 7,14) || herm(4,2) != cplx( 0,  0) || herm(4,3) != cplx(14,-2) || herm(4,4) != cplx( 1,  0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2,  0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,  0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (14,  6) (18,-9) (25,  0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) ( 7,-14) ( 8,-2) )\n"
                                        "( (14,-6) ( 0, 0) (18,  0) (11,11) ( 0,  0) (-2, 0) )\n"
                                        "( (18, 9) ( 0, 0) (11,-11) (20, 0) (14,  2) ( 0, 0) )\n"
                                        "( (25, 0) ( 7,14) ( 0,  0) (14,-2) ( 1,  0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4,  0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(-9, 0);
         mat(0,2) = cplx( 5, 0);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-7, 0);
         mat(2,0) = cplx( 6, 0);
         mat(2,1) = cplx(11, 0);
         mat(2,2) = cplx(99,99);
         mat(3,0) = cplx(11, 0);
         mat(3,1) = cplx( 4, 0);
         mat(3,2) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );
         sm %= mat;

         checkRows    ( herm,  6UL );
         checkColumns ( herm,  6UL );
         checkNonZeros( herm, 26UL );

         if( sm(0,0) != cplx(14,  6) || sm(0,1) != cplx(18,-9) || sm(0,2) != cplx(25,  0) ||
             sm(1,0) != cplx( 0,  0) || sm(1,1) != cplx( 0, 0) || sm(1,2) != cplx( 7,-14) ||
             sm(2,0) != cplx(18,  0) || sm(2,1) != cplx(11,11) || sm(2,2) != cplx( 0,  0) ||
             sm(3,0) != cplx(11,-11) || sm(3,1) != cplx(20, 0) || sm(3,2) != cplx(14,  2) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << sm << "\n"
                << "   Expected result:\n( (14,  6) (18,-9) (25,  0) )\n"
                                        "( ( 0,  0) ( 0, 0) ( 7,-14) )\n"
                                        "( (18,  0) (11,11) ( 0,  0) )\n"
                                        "( (11,-11) (20, 0) (14,  2) )\n";
            throw std::runtime_error( oss.str() );
         }

         if( herm(0,0) != cplx( 1, 0) || herm(0,1) != cplx(-4,-1) || herm(0,2) != cplx(14,  6) || herm(0,3) != cplx(18,-9) || herm(0,4) != cplx(25,  0) || herm(0,5) != cplx( 0, 0) ||
             herm(1,0) != cplx(-4, 1) || herm(1,1) != cplx( 2, 0) || herm(1,2) != cplx( 0,  0) || herm(1,3) != cplx( 0, 0) || herm(1,4) != cplx( 7,-14) || herm(1,5) != cplx( 8,-2) ||
             herm(2,0) != cplx(14,-6) || herm(2,1) != cplx( 0, 0) || herm(2,2) != cplx(18,  0) || herm(2,3) != cplx(11,11) || herm(2,4) != cplx( 0,  0) || herm(2,5) != cplx(-2, 0) ||
             herm(3,0) != cplx(18, 9) || herm(3,1) != cplx( 0, 0) || herm(3,2) != cplx(11,-11) || herm(3,3) != cplx(20, 0) || herm(3,4) != cplx(14,  2) || herm(3,5) != cplx( 0, 0) ||
             herm(4,0) != cplx(25, 0) || herm(4,1) != cplx( 7,14) || herm(4,2) != cplx( 0,  0) || herm(4,3) != cplx(14,-2) || herm(4,4) != cplx( 1,  0) || herm(4,5) != cplx(-4, 0) ||
             herm(5,0) != cplx( 0, 0) || herm(5,1) != cplx( 8, 2) || herm(5,2) != cplx(-2,  0) || herm(5,3) != cplx( 0, 0) || herm(5,4) != cplx(-4,  0) || herm(5,5) != cplx( 7, 0) ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment to submatrix failed\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n"
                << "   Expected result:\n( ( 1, 0) (-4,-1) (14,  6) (18,-9) (25,  0) ( 0, 0) )\n"
                                        "( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) ( 7,-14) ( 8,-2) )\n"
                                        "( (14,-6) ( 0, 0) (18,  0) (11,11) ( 0,  0) (-2, 0) )\n"
                                        "( (18, 9) ( 0, 0) (11,-11) (20, 0) (14,  2) ( 0, 0) )\n"
                                        "( (25, 0) ( 7,14) ( 0,  0) (14,-2) ( 1,  0) (-4, 0) )\n"
                                        "( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4,  0) ( 7, 0) )\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (11,  0) (20, 5) (28,12) (16,-8) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (24, -6) (12, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (28,-12) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (16,  8) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,  0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,  0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix Schur product assignment test 5";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-5, 0);
         mat(0,2) = cplx( 4, 0);
         mat(0,3) = cplx(-8, 0);
         mat(1,0) = cplx(-6, 0);
         mat(1,1) = cplx( 6, 0);
         mat(1,2) = cplx(99,99);
         mat(1,3) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-6, 0);
         mat(1,0) = cplx(-5, 0);
         mat(1,1) = cplx( 6, 0);
         mat(2,0) = cplx( 4, 0);
         mat(2,1) = cplx(99,99);
         mat(3,0) = cplx(-8, 0);
         mat(3,1) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,  3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) (18,  0) (14,14) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) (22,-22) (20, 0) (21, 3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0,  0) (21,-3) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix Schur product assignment test 6";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 6, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(99,99);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(22, 0);
         mat(1,2) = cplx( 4, 0);
         mat(1,3) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(99,99);
         mat(1,0) = cplx( 6, 0);
         mat(1,1) = cplx(22, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx( 4, 0);
         mat(3,0) = cplx(99,99);
         mat(3,1) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (16, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) (21,3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) (21,-3) (14,0) (20, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (16, 0) ( 0, 0) (24,0) (28, 0) )
   {
      test_ = "Sparse matrix Schur product assignment test 7";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 3, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(-5, 0);
         mat(1,0) = cplx(-8, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-6, 0);
         mat(1,3) = cplx( 4, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(-8, 0);
         mat(1,0) = cplx( 3, 0);
         mat(1,1) = cplx(99,99);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(-6, 0);
         mat(3,0) = cplx(-5, 0);
         mat(3,1) = cplx( 4, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (14,  6) (18,-9) (25, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) ( 7,14) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,-6) ( 0, 0) (18,  0) (11,11) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (18, 9) ( 0, 0) (22,-22) (20, 0) (14, 2) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (25, 0) ( 7,14) ( 0,  0) (14,-2) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix Schur product assignment test 8";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(99,99);
         mat(0,2) = cplx( 6, 0);
         mat(0,3) = cplx(11, 0);
         mat(1,0) = cplx(-9, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(22, 0);
         mat(1,3) = cplx( 4, 0);
         mat(2,0) = cplx( 5, 0);
         mat(2,1) = cplx(-7, 0);
         mat(2,2) = cplx(99,99);
         mat(2,3) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(-9, 0);
         mat(0,2) = cplx( 5, 0);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-7, 0);
         mat(2,0) = cplx( 6, 0);
         mat(2,1) = cplx(22, 0);
         mat(2,2) = cplx(99,99);
         mat(3,0) = cplx(11, 0);
         mat(3,1) = cplx( 4, 0);
         mat(3,2) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( (11,  0) (20, 5) (28,12) (16,-8) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (20, -5) (12, 4) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (28,-12) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (16,  8) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5,  0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0,  0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )
   {
      test_ = "Sparse matrix Schur product assignment test 9";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-5, 0);
         mat(0,2) = cplx( 4, 0);
         mat(0,3) = cplx(-8, 0);
         mat(1,0) = cplx(-5, 0);
         mat(1,1) = cplx( 6, 2);
         mat(1,2) = cplx(99,99);
         mat(1,3) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(11, 0);
         mat(0,1) = cplx(-5, 0);
         mat(1,0) = cplx(-5, 0);
         mat(1,1) = cplx( 6, 2);
         mat(2,0) = cplx( 4, 0);
         mat(2,1) = cplx(99,99);
         mat(3,0) = cplx(-8, 0);
         mat(3,1) = cplx(99,99);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 0UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7,  3) (-2, 1) ( 5, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) (-1, 2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) (18,  0) (14,14) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) (14,-14) (20, 5) (21, 3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0,  0) (21,-3) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix Schur product assignment test 10";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 6, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(99,99);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(14, 0);
         mat(1,2) = cplx( 4, 1);
         mat(1,3) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 1UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(99,99);
         mat(1,0) = cplx( 6, 0);
         mat(1,1) = cplx(14, 0);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx( 4, 1);
         mat(3,0) = cplx(99,99);
         mat(3,1) = cplx( 3, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 1UL, 2UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (16, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) (21,3) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( ( 5, 0) (-1,-2) ( 0, 0) (21,-3) (14,0) (20, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (16, 0) ( 0, 0) (20,0) (28, 7) )
   {
      test_ = "Sparse matrix Schur product assignment test 11";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 2UL, 4UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx( 3, 0);
         mat(0,2) = cplx(14, 0);
         mat(0,3) = cplx(-5, 0);
         mat(1,0) = cplx(-8, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-5, 0);
         mat(1,3) = cplx( 4, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 4UL, 2UL, 2UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 2UL, 8UL );
         mat(0,0) = cplx(99,99);
         mat(0,1) = cplx(-8, 0);
         mat(1,0) = cplx( 3, 0);
         mat(1,1) = cplx(99,99);
         mat(2,0) = cplx(14, 0);
         mat(2,1) = cplx(-5, 0);
         mat(3,0) = cplx(-5, 0);
         mat(3,1) = cplx( 4, 1);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 4UL, 4UL, 2UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }
   }

   // ( ( 1, 0) (-4,-1) ( 7, 3) (-2, 1) ( 5,0) ( 0, 0) )      ( ( 1, 0) (-4,-1) (14,  6) (18,-9) (25, 0) ( 0, 0) )
   // ( (-4, 1) ( 2, 0) ( 0, 0) ( 0, 0) (-1,2) ( 8,-2) )      ( (-4, 1) ( 2, 0) ( 0,  0) ( 0, 0) ( 7,14) ( 8,-2) )
   // ( ( 7,-3) ( 0, 0) ( 3, 0) ( 1, 1) ( 0,0) (-2, 0) )  =>  ( (14,-6) ( 0, 0) (18,  0) (11,11) ( 0, 0) (-2, 0) )
   // ( (-2,-1) ( 0, 0) ( 1,-1) ( 5, 0) ( 7,1) ( 0, 0) )      ( (18, 9) ( 0, 0) (11,-11) (20, 5) (14, 2) ( 0, 0) )
   // ( ( 5, 0) (-1,-2) ( 0, 0) ( 7,-1) ( 1,0) (-4, 0) )      ( (25, 0) ( 7,14) ( 0,  0) (14,-2) ( 1, 0) (-4, 0) )
   // ( ( 0, 0) ( 8, 2) (-2, 0) ( 0, 0) (-4,0) ( 7, 0) )      ( ( 0, 0) ( 8, 2) (-2,  0) ( 0, 0) (-4, 0) ( 7, 0) )
   {
      test_ = "Sparse matrix Schur product assignment test 12";

      {
         blaze::CompressedMatrix<cplx,blaze::rowMajor> mat( 3UL, 4UL, 12UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(99,99);
         mat(0,2) = cplx( 6, 0);
         mat(0,3) = cplx(11, 0);
         mat(1,0) = cplx(-9, 0);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(11, 0);
         mat(1,3) = cplx( 4, 1);
         mat(2,0) = cplx( 5, 0);
         mat(2,1) = cplx(-7, 0);
         mat(2,2) = cplx(99,99);
         mat(2,3) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 2UL, 0UL, 3UL, 4UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
            throw std::runtime_error( oss.str() );
         }
         catch( std::invalid_argument& ) {}
      }

      {
         blaze::CompressedMatrix<cplx,blaze::columnMajor> mat( 4UL, 3UL, 12UL );
         mat(0,0) = cplx( 2, 0);
         mat(0,1) = cplx(-9, 0);
         mat(0,2) = cplx( 5, 0);
         mat(1,0) = cplx(99,99);
         mat(1,1) = cplx(99,99);
         mat(1,2) = cplx(-7, 0);
         mat(2,0) = cplx( 6, 0);
         mat(2,1) = cplx(11, 0);
         mat(2,2) = cplx(99,99);
         mat(3,0) = cplx(11, 0);
         mat(3,1) = cplx( 4, 1);
         mat(3,2) = cplx( 2, 0);

         HT herm;
         init( herm );

         auto sm = submatrix( herm, 0UL, 2UL, 4UL, 3UL );

         try {
            sm %= mat;

            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment of invalid matrix succeeded\n"
                << " Details:\n"
                << "   Result:\n" << herm << "\n";
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
void SubmatrixComplexTest::checkRows( const Type& matrix, size_t expectedRows ) const
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
void SubmatrixComplexTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
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
void SubmatrixComplexTest::checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const
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
void SubmatrixComplexTest::init( HT& herm )
{
   herm.resize( 6UL );
   herm(0,0) = cplx( 1, 0);
   herm(0,1) = cplx(-4,-1);
   herm(0,2) = cplx( 7, 3);
   herm(0,3) = cplx(-2, 1);
   herm(0,4) = cplx( 5, 0);
   herm(1,1) = cplx( 2, 0);
   herm(1,4) = cplx(-1, 2);
   herm(1,5) = cplx( 8,-2);
   herm(2,2) = cplx( 3, 0);
   herm(2,3) = cplx( 1, 1);
   herm(2,5) = cplx(-2, 0);
   herm(3,3) = cplx( 5, 0);
   herm(3,4) = cplx( 7, 1);
   herm(4,4) = cplx( 1, 0);
   herm(4,5) = cplx(-4, 0);
   herm(5,5) = cplx( 7, 0);
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the assignment to a submatrix of a HermitianMatrix.
//
// \return void
*/
void runTest()
{
   SubmatrixComplexTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the HermitianMatrix submatrix complex test.
*/
#define RUN_HERMITIANMATRIX_SUBMATRIXCOMPLEX_TEST \
   blazetest::mathtest::adaptors::hermitianmatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace hermitianmatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
