//=================================================================================================
/*!
//  \file blazetest/mathtest/views/submatrix/DenseAlignedTest.h
//  \brief Header file for the Submatrix dense aligned test
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

#ifndef _BLAZETEST_MATHTEST_VIEWS_SUBMATRIX_DENSEALIGNEDTEST_H_
#define _BLAZETEST_MATHTEST_VIEWS_SUBMATRIX_DENSEALIGNEDTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Submatrix.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace views {

namespace submatrix {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the dense aligned Submatrix specialization.
//
// This class represents a test suite for the blaze::Submatrix class template specialization for
// dense aligned submatrices. It performs a series of both compile time as well as runtime tests.
*/
class DenseAlignedTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit DenseAlignedTest();
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
   void testConstructors();
   void testAssignment();
   void testAddAssign();
   void testSubAssign();
   void testSchurAssign();
   void testMultAssign();
   void testScaling();
   void testFunctionCall();
   void testIterator();
   void testNonZeros();
   void testReset();
   void testClear();
   void testTranspose();
   void testCTranspose();
   void testIsDefault();
   void testIsSame();
   void testSubmatrix();
   void testRow();
   void testRows();
   void testColumn();
   void testColumns();
   void testBand();

   template< typename Type >
   void checkRows( const Type& matrix, size_t expectedRows ) const;

   template< typename Type >
   void checkColumns( const Type& matrix, size_t expectedColumns ) const;

   template< typename Type >
   void checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const;

   template< typename Type >
   void checkNonZeros( const Type& matrix, size_t index, size_t expectedNonZeros ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initialize();
   //@}
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   using MT    = blaze::DynamicMatrix<int,blaze::rowMajor>;  //!< Row-major dynamic matrix type
   using OMT   = MT::OppositeType;                           //!< Column-major dynamic matrix type
   using ASMT  = blaze::Submatrix<MT,blaze::aligned>;        //!< Aligned dense submatrix type for row-major matrices.
   using USMT  = blaze::Submatrix<MT,blaze::unaligned>;      //!< Unaligned dense submatrix type for row-major matrices.
   using AOSMT = blaze::Submatrix<OMT,blaze::aligned>;       //!< Aligned dense submatrix type for column-major matrices.
   using UOSMT = blaze::Submatrix<OMT,blaze::unaligned>;     //!< Unaligned dense submatrix type for column-major matrices.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MT  mat1_;   //!< First row-major dynamic matrix.
                /*!< The \f$ 64 \times 64 \f$ row-major dense matrix is randomly initialized. */
   MT  mat2_;   //!< Second row-major dynamic matrix.
                /*!< The \f$ 64 \times 64 \f$ row-major dense matrix is randomly initialized. */
   OMT tmat1_;  //!< First column-major dynamic matrix.
                /*!< The \f$ 64 \times 64 \f$ column-major dense matrix is randomly initialized. */
   OMT tmat2_;  //!< Second column-major dynamic matrix.
                /*!< The \f$ 64 \times 64 \f$ column-major dense matrix is randomly initialized. */

   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OMT   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ASMT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( USMT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( AOSMT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( UOSMT );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checking the number of rows of the given dense matrix.
//
// \param matrix The dense matrix to be checked.
// \param expectedRows The expected number of rows of the dense matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of rows of the given dense matrix. In case the
// actual number of rows does not correspond to the given expected number of rows, a
// \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense matrix
void DenseAlignedTest::checkRows( const Type& matrix, size_t expectedRows ) const
{
   if( rows( matrix ) != expectedRows ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of rows detected\n"
          << " Details:\n"
          << "   Number of rows         : " << rows( matrix ) << "\n"
          << "   Expected number of rows: " << expectedRows << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of columns of the given dense matrix.
//
// \param matrix The dense matrix to be checked.
// \param expectedColumns The expected number of columns of the dense matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of columns of the given dense matrix. In case the
// actual number of columns does not correspond to the given expected number of columns,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense matrix
void DenseAlignedTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
{
   if( columns( matrix ) != expectedColumns ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of columns detected\n"
          << " Details:\n"
          << "   Number of columns         : " << columns( matrix ) << "\n"
          << "   Expected number of columns: " << expectedColumns << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of non-zero elements of the given dense matrix.
//
// \param matrix The dense matrix to be checked.
// \param expectedNonZeros The expected number of non-zero elements of the dense matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given dense matrix. In
// case the actual number of non-zero elements does not correspond to the given expected
// number, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense matrix
void DenseAlignedTest::checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const
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


//*************************************************************************************************
/*!\brief Checking the number of non-zero elements in a specific row/column of the given dense matrix.
//
// \param matrix The dense matrix to be checked.
// \param index The row/column to be checked.
// \param expectedNonZeros The expected number of non-zero elements in the specified row/column.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements in the specified row/column of the
// given dense matrix. In case the actual number of non-zero elements does not correspond
// to the given expected number, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense matrix
void DenseAlignedTest::checkNonZeros( const Type& matrix, size_t index, size_t expectedNonZeros ) const
{
   if( nonZeros( matrix, index ) != expectedNonZeros ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of non-zero elements in "
          << ( blaze::IsRowMajorMatrix<Type>::value ? "row " : "column " ) << index << "\n"
          << " Details:\n"
          << "   Number of non-zeros         : " << nonZeros( matrix, index ) << "\n"
          << "   Expected number of non-zeros: " << expectedNonZeros << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( capacity( matrix, index ) < nonZeros( matrix, index ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected in "
          << ( blaze::IsRowMajorMatrix<Type>::value ? "row " : "column " ) << index << "\n"
          << " Details:\n"
          << "   Number of non-zeros: " << nonZeros( matrix, index ) << "\n"
          << "   Capacity           : " << capacity( matrix, index ) << "\n";
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
/*!\brief Testing the functionality of the dense aligned Submatrix specialization.
//
// \return void
*/
void runTest()
{
   DenseAlignedTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the Submatrix dense aligned test.
*/
#define RUN_SUBMATRIX_DENSEALIGNED_TEST \
   blazetest::mathtest::views::submatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace submatrix

} // namespace views

} // namespace mathtest

} // namespace blazetest

#endif
