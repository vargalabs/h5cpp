//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/diagonalmatrix/SparseTest.h
//  \brief Header file for the DiagonalMatrix sparse test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_DIAGONALMATRIX_SPARSETEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_DIAGONALMATRIX_SPARSETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Diagonal.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace diagonalmatrix {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the sparse DiagonalMatrix specialization.
//
// This class represents a test suite for the blaze::DiagonalMatrix class template specialization
// for sparse matrices. It performs a series of both compile time as well as runtime tests.
*/
class SparseTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit SparseTest();
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
   void testAssignment  ();
   void testAddAssign   ();
   void testSubAssign   ();
   void testSchurAssign ();
   void testMultAssign  ();
   void testScaling     ();
   void testFunctionCall();
   void testIterator    ();
   void testNonZeros    ();
   void testReset       ();
   void testClear       ();
   void testResize      ();
   void testReserve     ();
   void testTrim        ();
   void testShrinkToFit ();
   void testSwap        ();
   void testSet         ();
   void testInsert      ();
   void testAppend      ();
   void testErase       ();
   void testFind        ();
   void testLowerBound  ();
   void testUpperBound  ();
   void testIsDefault   ();
   void testSubmatrix   ();
   void testRow         ();
   void testColumn      ();

   template< typename Type >
   void checkRows( const Type& matrix, size_t expectedRows ) const;

   template< typename Type >
   void checkColumns( const Type& matrix, size_t expectedColumns ) const;

   template< typename Type >
   void checkCapacity( const Type& matrix, size_t minCapacity ) const;

   template< typename Type >
   void checkCapacity( const Type& matrix, size_t index, size_t minCapacity ) const;

   template< typename Type >
   void checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const;

   template< typename Type >
   void checkNonZeros( const Type& matrix, size_t index, size_t expectedNonZeros ) const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Type of the row-major diagonal matrix.
   using DT = blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> >;

   //! Type of the column-major diagonal matrix.
   using ODT = blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> >;

   using RDT  = DT::Rebind<double>::Other;   //!< Rebound row-major diagonal matrix type.
   using ORDT = ODT::Rebind<double>::Other;  //!< Rebound column-major diagonal matrix type.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( DT                  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( DT::ResultType      );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( DT::OppositeType    );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( DT::TransposeType   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ODT                 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ODT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ODT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ODT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RDT                 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RDT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RDT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RDT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ORDT                );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ORDT::ResultType    );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ORDT::OppositeType  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ORDT::TransposeType );

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( DT                  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( DT::ResultType      );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( DT::OppositeType    );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( DT::TransposeType   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ODT                 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ODT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( ODT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( ODT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RDT                 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RDT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( RDT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( RDT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ORDT                );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ORDT::ResultType    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( ORDT::OppositeType  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( ORDT::TransposeType );

   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( DT                  );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( DT::ResultType      );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( DT::OppositeType    );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( DT::TransposeType   );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( ODT                 );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( ODT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( ODT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( ODT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( RDT                 );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( RDT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( RDT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( RDT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( ORDT                );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( ORDT::ResultType    );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( ORDT::OppositeType  );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE( ORDT::TransposeType );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DT::ResultType      );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DT::OppositeType    );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DT::TransposeType   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ODT::ResultType     );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ODT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ODT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RDT::ResultType     );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RDT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RDT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ORDT::ResultType    );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ORDT::OppositeType  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ORDT::TransposeType );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DT::ElementType,   DT::ResultType::ElementType      );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DT::ElementType,   DT::OppositeType::ElementType    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DT::ElementType,   DT::TransposeType::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ODT::ElementType,  ODT::ResultType::ElementType     );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ODT::ElementType,  ODT::OppositeType::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ODT::ElementType,  ODT::TransposeType::ElementType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RDT::ElementType,  RDT::ResultType::ElementType     );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RDT::ElementType,  RDT::OppositeType::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RDT::ElementType,  RDT::TransposeType::ElementType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ORDT::ElementType, ORDT::ResultType::ElementType    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ORDT::ElementType, ORDT::OppositeType::ElementType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ORDT::ElementType, ORDT::TransposeType::ElementType );
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
void SparseTest::checkRows( const Type& matrix, size_t expectedRows ) const
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
// \param expectedRows The expected number of columns of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of columns of the given matrix. In case the actual number of
// columns does not correspond to the given expected number of columns, a \a std::runtime_error
// exception is thrown.
*/
template< typename Type >  // Type of the matrix
void SparseTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
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
/*!\brief Checking the capacity of the given matrix.
//
// \param matrix The matrix to be checked.
// \param minCapacity The expected minimum capacity of the matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given matrix. In case the actual capacity is smaller
// than the given expected minimum capacity, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the matrix
void SparseTest::checkCapacity( const Type& matrix, size_t minCapacity ) const
{
   if( capacity( matrix ) < minCapacity ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Capacity                 : " << capacity( matrix ) << "\n"
          << "   Expected minimum capacity: " << minCapacity << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the capacity of a specific row/column of the given matrix.
//
// \param matrix The matrix to be checked.
// \param index The row/column to be checked.
// \param minCapacity The expected minimum capacity of the specified row/column.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of a specific row/column of the given matrix. In case the
// actual capacity is smaller than the given expected minimum capacity, a \a std::runtime_error
// exception is thrown.
*/
template< typename Type >  // Type of the matrix
void SparseTest::checkCapacity( const Type& matrix, size_t index, size_t minCapacity ) const
{
   if( capacity( matrix, index ) < minCapacity ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected in "
          << ( blaze::IsRowMajorMatrix<Type>::value ? "row " : "column " ) << index << "\n"
          << " Details:\n"
          << "   Capacity                 : " << capacity( matrix, index ) << "\n"
          << "   Expected minimum capacity: " << minCapacity << "\n";
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
void SparseTest::checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const
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
/*!\brief Checking the number of non-zero elements in a specific row/column of the given matrix.
//
// \param matrix The matrix to be checked.
// \param index The row/column to be checked.
// \param expectedNonZeros The expected number of non-zero elements in the specified row/column.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements in the specified row/column of the
// given matrix. In case the actual number of non-zero elements does not correspond to the
// given expected number, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the matrix
void SparseTest::checkNonZeros( const Type& matrix, size_t index, size_t expectedNonZeros ) const
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
/*!\brief Testing the functionality of the sparse DiagonalMatrix specialization.
//
// \return void
*/
void runTest()
{
   SparseTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the DiagonalMatrix sparse test.
*/
#define RUN_DIAGONALMATRIX_SPARSE_TEST \
   blazetest::mathtest::adaptors::diagonalmatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace diagonalmatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
