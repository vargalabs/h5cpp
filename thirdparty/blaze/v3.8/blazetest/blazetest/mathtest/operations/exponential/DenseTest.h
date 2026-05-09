//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/exponential/DenseTest.h
//  \brief Header file for the dense matrix exponential test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_EXPONENTIAL_DENSETEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_EXPONENTIAL_DENSETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/Aliases.h>
#include <blaze/math/DenseMatrix.h>
#include <blaze/math/Submatrix.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/IsEqual.h>
#include <blazetest/system/LAPACK.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace operations {

namespace exponential {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all dense matrix exponential tests.
//
// This class represents a test suite for the dense matrix exponential functionality. It performs
// a series of matrix exponentials on all dense matrix types of the Blaze library.
*/
class DenseTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit DenseTest();
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
   void testSpecific();

   template< typename Type >
   void testRandom( size_t N );
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
/*!\brief Test of the exponential functionality with random \f$ N \times N \f$ matrices.
//
// \param N The number of rows and columns of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix exponential for random \f$ N \times N \f$ matrices. In
// case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom( size_t N )
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using namespace blaze;

   {
      test_ = "Matrix exponential";

      Type A;
      resize( A, N, N );
      randomize( A );

      OppositeType_t<Type> B( A );

      const auto ExpA( evaluate( matexp( A ) ) );
      const auto ExpB( evaluate( matexp( B ) ) );

      if( !isEqual( ExpA, ExpB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Matrix exponential failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ElementType_t<Type> ).name() << "\n"
             << "   Initial row-major matrix:\n" << A << "\n"
             << "   Initial column-major matrix:\n" << B << "\n"
             << "   Row-major matrix exponential:\n" << ExpA << "\n"
             << "   Column-major matrix exponential:\n" << ExpB << "\n";
         throw std::runtime_error( oss.str() );
      }

      const auto ExpTA( evaluate( matexp( trans( A ) ) ) );
      const auto ExpTB( evaluate( matexp( trans( B ) ) ) );

      if( !isEqual( ExpTA, ExpTB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose matrix exponential failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ElementType_t<Type> ).name() << "\n"
             << "   Initial row-major matrix:\n" << trans( B ) << "\n"
             << "   Initial column-major matrix:\n" << trans( A ) << "\n"
             << "   Row-major matrix exponential:\n" << ExpTB << "\n"
             << "   Column-major matrix exponential:\n" << ExpTA << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Submatrix exponential";

      Type A;
      resize( A, N, N );
      randomize( A );

      OppositeType_t<Type> B( A );

      const auto ExpA( evaluate( matexp( submatrix( A, 0UL, 0UL, N, N ) ) ) );
      const auto ExpB( evaluate( matexp( submatrix( B, 0UL, 0UL, N, N ) ) ) );

      if( !isEqual( ExpA, ExpB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Submatrix exponential failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ElementType_t<Type> ).name() << "\n"
             << "   Initial row-major matrix:\n" << A << "\n"
             << "   Initial column-major matrix:\n" << B << "\n"
             << "   Row-major matrix exponential:\n" << ExpA << "\n"
             << "   Column-major matrix exponential:\n" << ExpB << "\n";
         throw std::runtime_error( oss.str() );
      }

      const auto ExpTA( evaluate( matexp( trans( submatrix( A, 0UL, 0UL, N, N ) ) ) ) );
      const auto ExpTB( evaluate( matexp( trans( submatrix( B, 0UL, 0UL, N, N ) ) ) ) );

      if( !isEqual( ExpTA, ExpTB ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Transpose submatrix exponential failed\n"
             << " Details:\n"
             << "   Matrix type:\n"
             << "     " << typeid( Type ).name() << "\n"
             << "   Element type:\n"
             << "     " << typeid( ElementType_t<Type> ).name() << "\n"
             << "   Initial row-major matrix:\n" << trans( B ) << "\n"
             << "   Initial column-major matrix:\n" << trans( A ) << "\n"
             << "   Row-major matrix exponential:\n" << ExpTB << "\n"
             << "   Column-major matrix exponential:\n" << ExpTA << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the dense matrix exponential.
//
// \return void
*/
void runTest()
{
   DenseTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the dense matrix exponential test.
*/
#define RUN_EXPONENTIAL_DENSE_TEST \
   blazetest::mathtest::operations::exponential::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace exponential

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
