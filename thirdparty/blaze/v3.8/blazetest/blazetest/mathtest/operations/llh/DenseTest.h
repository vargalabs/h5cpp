//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/llh/DenseTest.h
//  \brief Header file for the dense matrix LLH test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_LLH_DENSETEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_LLH_DENSETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/Aliases.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Random.h>
#include <blazetest/system/LAPACK.h>


namespace blazetest {

namespace mathtest {

namespace operations {

namespace llh {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all dense matrix LLH tests.
//
// This class represents a test suite for the dense matrix LLH decomposition functionality. It
// performs a series of LLH decompositions on all dense matrix types of the Blaze library.
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
   template< typename Type >
   void testRandom();

   void testGeneral();
   void testSymmetric();
   void testHermitian();
   void testLower();
   void testUniLower();
   void testUpper();
   void testUniUpper();
   void testDiagonal();
   //@}
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   using cfloat  = blaze::complex<float>;   //!< Single precision complex test type.
   using cdouble = blaze::complex<double>;  //!< Double precision complex test type.
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
/*!\brief Test of the LLH decomposition with a randomly initialized matrix of the given type.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function tests the dense matrix LLH decomposition for a randomly initialized matrix of the
// given type. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   test_ = "LLH decomposition";

   using MT = blaze::RemoveAdaptor_t<Type>;

   const size_t n( blaze::rand<size_t>( 3UL, 8UL ) );

   Type A;
   blaze::LowerMatrix<MT> L;

   resize( A, n, n );
   makePositiveDefinite( A );
   blaze::llh( A, L );

   const MT LLH( L * ctrans( L ) );

   if( LLH != A ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: LLH decomposition failed\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( blaze::ElementType_t<Type> ).name() << "\n"
          << "   Result:\n" << LLH << "\n"
          << "   Expected result:\n" << A << "\n";
      throw std::runtime_error( oss.str() );
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
/*!\brief Testing the dense matrix LLH decomposition.
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
/*!\brief Macro for the execution of the dense matrix LLH test.
*/
#define RUN_DENSE_LLH_TEST \
   blazetest::mathtest::operations::llh::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace llh

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
