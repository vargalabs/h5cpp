//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/eigen/DenseTest.h
//  \brief Header file for the dense matrix eigenvalue test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_EIGEN_DENSETEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_EIGEN_DENSETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/util/Random.h>


namespace blazetest {

namespace mathtest {

namespace operations {

namespace eigen {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all dense matrix eigenvalue/eigenvector tests.
//
// This class represents a test suite for the dense matrix eigenvalue/eigenvector functionality.
// It performs a series of eigenvalue/eigenvector computations on several dense matrix types of
// the Blaze library.
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
   void testGeneral();
   void testSymmetric();
   void testHermitian();
   void testLower();
   void testUpper();
   void testDiagonal();
   //@}
   //**********************************************************************************************

   //**Error detection functions*******************************************************************
   /*!\name Error detection functions */
   //@{
   template< typename VT, typename MT, bool SO, typename ST >
   void checkEigenvector( const blaze::DenseVector<VT,false>& v,
                          const blaze::DenseMatrix<MT,SO>& A, ST w );

   template< typename VT, typename MT, bool SO, typename ST >
   void checkEigenvector( const blaze::DenseVector<VT,true>& u,
                          const blaze::DenseMatrix<MT,SO>& A, ST w );
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
//  ERROR DETECTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checking the given right eigenvector.
//
// \param v The right eigenvector to be checked.
// \param A The corresponding dense matrix.
// \param w The corresponding eigenvalue.
// \return void
// \exception std::runtime_error Invalid left eigenvector detected.
//
// This function checks the given right eigenvector \f$v[j]\f$ by testing if it satisfies

                          \f[ A * v[j] = lambda[j] * v[j], \f]

// where \f$lambda[j]\f$ is the corresponding eigenvalue.
*/
template< typename VT    // Type of the eigenvector v
        , typename MT    // Type of the matrix A
        , bool SO        // Storage order of the matrix A
        , typename ST >  // Type of the eigenvalue w
void DenseTest::checkEigenvector( const blaze::DenseVector<VT,false>& v,
                                  const blaze::DenseMatrix<MT,SO>& A, ST w )
{
   if( A * v != w * v ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid right eigenvector detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   System matrix:\n" << A << "\n"
          << "   Eigenvalue = " << w << "\n"
          << "   Right eigenvector:\n" << v << "\n"
          << "   A * v =\n" << ( A * v ) << "\n"
          << "   w * v =\n" << ( w * v ) << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the given left eigenvector.
//
// \param u The left eigenvector to be checked.
// \param A The corresponding dense matrix.
// \param w The corresponding eigenvalue.
// \return void
// \exception std::runtime_error Invalid left eigenvector detected.
//
// This function checks the given left eigenvector \f$u[j]\f$ by testing if it satisfies

                       \f[ u[j]^{H} * A = lambda[j] * u[j]^{H}, \f]

// where \f$lambda[j]\f$ is the corresponding eigenvalue.
*/
template< typename VT    // Type of the eigenvector u
        , typename MT    // Type of the matrix A
        , bool SO        // Storage order of the matrix A
        , typename ST >  // Type of the eigenvalue w
void DenseTest::checkEigenvector( const blaze::DenseVector<VT,true>& u,
                                  const blaze::DenseMatrix<MT,SO>& A, ST w )
{
   if( u * A != u * w ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid left eigenvector detected\n"
          << " Details:\n"
          << "   Random seed = " << blaze::getSeed() << "\n"
          << "   System matrix:\n" << A << "\n"
          << "   Eigenvalue = " << w << "\n"
          << "   Left eigenvector:\n" << u << "\n"
          << "   v * A =\n" << ( u * A ) << "\n"
          << "   v * w =\n" << ( u * w ) << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the dense matrix eigenvalue functionality.
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
/*!\brief Macro for the execution of the dense matrix eigenvalue test.
*/
#define RUN_DENSE_EIGEN_TEST \
   blazetest::mathtest::operations::eigen::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace eigen

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
