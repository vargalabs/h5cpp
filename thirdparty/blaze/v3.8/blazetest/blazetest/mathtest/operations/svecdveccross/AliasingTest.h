//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/svecdveccross/AliasingTest.h
//  \brief Header file for the sparse vector/dense vector cross product aliasing test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_SVECDVECCROSS_ALIASINGTEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_SVECDVECCROSS_ALIASINGTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/StaticVector.h>


namespace blazetest {

namespace mathtest {

namespace operations {

namespace svecdveccross {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the sparse vector/dense vector cross product aliasing test.
//
// This class represents a test suite for all sparse vector/dense vector cross product aliasing
// tests. It performs a series of runtime tests to assure that all mathematical operations work
// correctly even in the presence of aliasing.
*/
class AliasingTest
{
 private:
   //**Type definitions****************************************************************************
   using DVec = blaze::DynamicVector<int,blaze::columnVector>;     //!< Dense column vector type.
   using DMat = blaze::DynamicMatrix<int,blaze::rowMajor>;         //!< Row-major dense matrix type.
   using SVec = blaze::CompressedVector<int,blaze::columnVector>;  //!< Sparse column vector type.
   using SMat = blaze::CompressedMatrix<int,blaze::rowMajor>;      //!< Row-major sparse matrix type.
   using RVec = blaze::StaticVector<int,3UL,blaze::columnVector>;  //!< Result column vector type.
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit AliasingTest();
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
   void testSVecDVecCross();

   template< typename T1, typename T2 >
   void checkResult( const T1& computedResult, const T2& expectedResult );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initialize();
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   SVec sa4_;     //!< The first sparse column vector.
                  /*!< The 4-dimensional vector is initialized as
                       \f[\left(\begin{array}{*{1}{c}}
                       -1 \\
                        0 \\
                       -3 \\
                        2 \\
                       \end{array}\right)\f]. */
   SVec sb3_;     //!< The second sparse column vector.
                  /*!< The 3-dimensional vector is initialized as
                       \f[\left(\begin{array}{*{1}{c}}
                       1 \\
                       2 \\
                       3 \\
                       \end{array}\right)\f]. */
   DVec da4_;     //!< The first dense column vector.
                  /*!< The 4-dimensional vector is initialized as
                       \f[\left(\begin{array}{*{1}{c}}
                       -1 \\
                        0 \\
                       -3 \\
                        2 \\
                       \end{array}\right)\f]. */
   DVec db3_;     //!< The second dense column vector.
                  /*!< The 3-dimensional vector is initialized as
                       \f[\left(\begin{array}{*{1}{c}}
                       1 \\
                       2 \\
                       3 \\
                       \end{array}\right)\f]. */
   SMat sA3x4_;   //!< The first row-major sparse matrix.
                  /*!< The \f$ 3 \times 4 \f$ matrix is initialized as
                       \f[\left(\begin{array}{*{3}{c}}
                       -1 & 0 & -2 & 0 \\
                        0 & 2 & -3 & 1 \\
                        0 & 1 &  2 & 2 \\
                       \end{array}\right)\f]. */
   SMat sB3x3_;   //!< The second row-major sparse matrix.
                  /*!< The \f$ 3 \times 3 \f$ matrix is initialized as
                       \f[\left(\begin{array}{*{3}{c}}
                       0 & -1 &  0 \\
                       1 & -2 &  2 \\
                       0 &  0 & -3 \\
                       \end{array}\right)\f]. */
   DMat dA3x4_;   //!< The first row-major dense matrix.
                  /*!< The \f$ 3 \times 4 \f$ matrix is initialized as
                       \f[\left(\begin{array}{*{3}{c}}
                       -1 & 0 & -2 & 0 \\
                        0 & 2 & -3 & 1 \\
                        0 & 1 &  2 & 2 \\
                       \end{array}\right)\f]. */
   DMat dB3x3_;   //!< The second row-major dense matrix.
                  /*!< The \f$ 3 \times 3 \f$ matrix is initialized as
                       \f[\left(\begin{array}{*{3}{c}}
                       0 & -1 &  0 \\
                       1 & -2 &  2 \\
                       0 &  0 & -3 \\
                       \end{array}\right)\f]. */
   RVec result_;  //!< The dense vector for the reference result.

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
/*!\brief Checking and comparing the computed result.
//
// \param computedResult The computed result.
// \param expectedResult The expected result.
// \return void
// \exception std::runtime_error Incorrect result detected.
//
// This function is called after each test case to check and compare the computed result.
// In case the computed and the expected result differ in any way, a \a std::runtime_error
// exception is thrown.
*/
template< typename T1    // Vector type of the computed result
        , typename T2 >  // Vector type of the expected result
void AliasingTest::checkResult( const T1& computedResult, const T2& expectedResult )
{
   if( computedResult != expectedResult ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Incorrect result detected\n"
          << " Details:\n"
          << "   Computed result:\n" << computedResult << "\n"
          << "   Expected result:\n" << expectedResult << "\n";
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
/*!\brief Testing the sparse vector/dense vector cross product in the presence of aliasing.
//
// \return void
*/
void runTest()
{
   AliasingTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the sparse vector/dense vector cross product aliasing test.
*/
#define RUN_SVECDVECCROSS_ALIASING_TEST \
   blazetest::mathtest::operations::svecdveccross::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace svecdveccross

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
