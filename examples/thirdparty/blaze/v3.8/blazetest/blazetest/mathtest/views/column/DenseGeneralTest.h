//=================================================================================================
/*!
//  \file blazetest/mathtest/views/column/DenseGeneralTest.h
//  \brief Header file for the Column dense general test
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

#ifndef _BLAZETEST_MATHTEST_VIEWS_COLUMN_DENSEGENERALTEST_H_
#define _BLAZETEST_MATHTEST_VIEWS_COLUMN_DENSEGENERALTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Column.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace views {

namespace column {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the dense general Column specialization.
//
// This class represents a test suite for the blaze::Column class template specialization for
// dense general matrices. It performs a series of both compile time as well as runtime tests.
*/
class DenseGeneralTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit DenseGeneralTest();
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
   void testMultAssign();
   void testDivAssign();
   void testCrossAssign();
   void testScaling();
   void testSubscript();
   void testIterator();
   void testNonZeros();
   void testReset();
   void testClear();
   void testIsDefault();
   void testIsSame();
   void testSubvector();
   void testElements();

   template< typename Type >
   void checkSize( const Type& column, size_t expectedSize ) const;

   template< typename Type >
   void checkRows( const Type& matrix, size_t expectedRows ) const;

   template< typename Type >
   void checkColumns( const Type& matrix, size_t expectedColumns ) const;

   template< typename Type >
   void checkCapacity( const Type& object, size_t minCapacity ) const;

   template< typename Type >
   void checkNonZeros( const Type& object, size_t expectedNonZeros ) const;

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
   using MT  = blaze::DynamicMatrix<int,blaze::rowMajor>;  //!< Row-major dynamic matrix type
   using OMT = MT::OppositeType;                           //!< Column-major dynamic matrix type
   using CT  = blaze::Column<MT>;                          //!< Dense column type for row-major matrices.
   using OCT = blaze::Column<OMT>;                         //!< Dense column type for column-major matrices.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MT  mat_;   //!< Row-major dynamic matrix.
               /*!< The \f$ 4 \times 5 \f$ matrix is initialized as
                    \f[\left(\begin{array}{*{5}{c}}
                    0 & 0 & -2 &  0 &  7 \\
                    0 & 1 &  0 &  4 & -8 \\
                    0 & 0 & -3 &  5 &  9 \\
                    0 & 0 &  0 & -6 & 10 \\
                    \end{array}\right)\f]. */
   OMT tmat_;  //!< Column-major dynamic matrix.
               /*!< The \f$ 4 \times 5 \f$ matrix is initialized as
                    \f[\left(\begin{array}{*{5}{c}}
                    0 & 0 & -2 &  0 &  7 \\
                    0 & 1 &  0 &  4 & -8 \\
                    0 & 0 & -3 &  5 &  9 \\
                    0 & 0 &  0 & -6 & 10 \\
                    \end{array}\right)\f]. */

   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( OMT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( CT  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( OCT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( CT  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( OCT );
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
/*!\brief Checking the size of the given dense column.
//
// \param column The dense column to be checked.
// \param expectedSize The expected size of the dense column.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the size of the given dense column. In case the actual size does not
// correspond to the given expected size, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense column
void DenseGeneralTest::checkSize( const Type& column, size_t expectedSize ) const
{
   if( size( column ) != expectedSize ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid size detected\n"
          << " Details:\n"
          << "   Size         : " << size( column ) << "\n"
          << "   Expected size: " << expectedSize << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of rows of the given dynamic matrix.
//
// \param matrix The dynamic matrix to be checked.
// \param expectedRows The expected number of rows of the dynamic matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of rows of the given dynamic matrix. In case the actual number
// of rows does not correspond to the given expected number of rows, a \a std::runtime_error
// exception is thrown.
*/
template< typename Type >  // Type of the dynamic matrix
void DenseGeneralTest::checkRows( const Type& matrix, size_t expectedRows ) const
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
/*!\brief Checking the number of columns of the given dynamic matrix.
//
// \param matrix The dynamic matrix to be checked.
// \param expectedColumns The expected number of columns of the dynamic matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of columns of the given dynamic matrix. In case the
// actual number of columns does not correspond to the given expected number of columns,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dynamic matrix
void DenseGeneralTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
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
/*!\brief Checking the capacity of the given dense column or dynamic matrix.
//
// \param object The dense column or dynamic matrix to be checked.
// \param minCapacity The expected minimum capacity.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given dense column or dynamic matrix. In
// case the actual capacity is smaller than the given expected minimum capacity, a
// \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense row or dynamic matrix
void DenseGeneralTest::checkCapacity( const Type& object, size_t minCapacity ) const
{
   if( capacity( object ) < minCapacity ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Capacity                 : " << capacity( object ) << "\n"
          << "   Expected minimum capacity: " << minCapacity << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of non-zero elements of the given dense column or dynamic matrix.
//
// \param object The dense column or dynamic matrix to be checked.
// \param expectedNonZeros The expected number of non-zero elements.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given dense column. In case
// the actual number of non-zero elements does not correspond to the given expected number,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dense column or dynamic matrix
void DenseGeneralTest::checkNonZeros( const Type& object, size_t expectedNonZeros ) const
{
   if( nonZeros( object ) != expectedNonZeros ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of non-zero elements\n"
          << " Details:\n"
          << "   Number of non-zeros         : " << nonZeros( object ) << "\n"
          << "   Expected number of non-zeros: " << expectedNonZeros << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( capacity( object ) < nonZeros( object ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Number of non-zeros: " << nonZeros( object ) << "\n"
          << "   Capacity           : " << capacity( object ) << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of non-zero elements in a specific row/column of the given dynamic matrix.
//
// \param matrix The dynamic matrix to be checked.
// \param index The row/column to be checked.
// \param expectedNonZeros The expected number of non-zero elements in the specified row/column.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements in the specified row/column of the
// given dynamic matrix. In case the actual number of non-zero elements does not correspond
// to the given expected number, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the dynamic matrix
void DenseGeneralTest::checkNonZeros( const Type& matrix, size_t index, size_t expectedNonZeros ) const
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
/*!\brief Testing the functionality of the dense general Column specialization.
//
// \return void
*/
void runTest()
{
   DenseGeneralTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the Column dense general test.
*/
#define RUN_COLUMN_DENSEGENERAL_TEST \
   blazetest::mathtest::views::column::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace column

} // namespace views

} // namespace mathtest

} // namespace blazetest

#endif
