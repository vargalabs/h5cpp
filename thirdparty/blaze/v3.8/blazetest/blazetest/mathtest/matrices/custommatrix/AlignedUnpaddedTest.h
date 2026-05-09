//=================================================================================================
/*!
//  \file blazetest/mathtest/matrices/custommatrix/AlignedUnpaddedTest.h
//  \brief Header file for the aligned/unpadded CustomMatrix class test
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

#ifndef _BLAZETEST_MATHTEST_MATRICES_CUSTOMMATRIX_ALIGNEDUNPADDEDTEST_H_
#define _BLAZETEST_MATHTEST_MATRICES_CUSTOMMATRIX_ALIGNEDUNPADDEDTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/CustomMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace custommatrix {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the aligned/unpadded CustomMatrix class template.
//
// This class represents a test suite for the specialization of the blaze::CustomMatrix class
// template for aligned and unpadded custom matrices. It performs a series of both compile
// time as well as runtime tests.
*/
class AlignedUnpaddedTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit AlignedUnpaddedTest();
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
   void testAt          ();
   void testIterator    ();
   void testNonZeros    ();
   void testReset       ();
   void testClear       ();
   void testSwap        ();
   void testTranspose   ();
   void testCTranspose  ();
   void testIsDefault   ();

   template< typename Type >
   void checkRows( const Type& matrix, size_t expectedRows ) const;

   template< typename Type >
   void checkColumns( const Type& matrix, size_t expectedColumns ) const;

   template< typename Type >
   void checkCapacity( const Type& matrix, size_t minCapacity ) const;

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
   //! Type of the row-major custom matrix.
   using MT = blaze::CustomMatrix<int,blaze::aligned,blaze::unpadded,blaze::rowMajor>;

   //! Type of the column-major custom matrix.
   using OMT = blaze::CustomMatrix<int,blaze::aligned,blaze::unpadded,blaze::columnMajor>;

   using RMT  = MT::Rebind<const double>::Other;   //!< Rebound row-major custom matrix type.
   using ORMT = OMT::Rebind<const double>::Other;  //!< Rebound column-major custom matrix type.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT                  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT::ResultType      );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT::OppositeType    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT::TransposeType   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OMT                 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OMT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OMT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OMT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RMT                 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RMT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RMT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RMT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ORMT                );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ORMT::ResultType    );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ORMT::OppositeType  );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ORMT::TransposeType );

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MT                  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( MT::ResultType      );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT::OppositeType    );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT::TransposeType   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OMT                 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OMT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( OMT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( OMT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RMT                 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( RMT::ResultType     );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( RMT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( RMT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ORMT                );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ORMT::ResultType    );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( ORMT::OppositeType  );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE   ( ORMT::TransposeType );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT::ResultType      );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT::OppositeType    );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT::TransposeType   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OMT::ResultType     );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OMT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OMT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RMT::ResultType     );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RMT::OppositeType   );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RMT::TransposeType  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ORMT::ResultType    );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ORMT::OppositeType  );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ORMT::TransposeType );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT::ElementType,   MT::ResultType::ElementType      );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT::ElementType,   MT::OppositeType::ElementType    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MT::ElementType,   MT::TransposeType::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( OMT::ElementType,  OMT::ResultType::ElementType     );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( OMT::ElementType,  OMT::OppositeType::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( OMT::ElementType,  OMT::TransposeType::ElementType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RMT::ElementType,  RMT::ResultType::ElementType     );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RMT::ElementType,  RMT::OppositeType::ElementType   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RMT::ElementType,  RMT::TransposeType::ElementType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ORMT::ElementType, ORMT::ResultType::ElementType    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ORMT::ElementType, ORMT::OppositeType::ElementType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ORMT::ElementType, ORMT::TransposeType::ElementType );
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
/*!\brief Checking the number of rows of the given custom matrix.
//
// \param matrix The custom matrix to be checked.
// \param expectedRows The expected number of rows of the custom matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of rows of the given custom matrix. In case the actual number
// of rows does not correspond to the given expected number of rows, a \a std::runtime_error
// exception is thrown.
*/
template< typename Type >  // Type of the custom matrix
void AlignedUnpaddedTest::checkRows( const Type& matrix, size_t expectedRows ) const
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
/*!\brief Checking the number of columns of the given custom matrix.
//
// \param matrix The custom matrix to be checked.
// \param expectedRows The expected number of columns of the custom matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of columns of the given custom matrix. In case the
// actual number of columns does not correspond to the given expected number of columns,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the custom matrix
void AlignedUnpaddedTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
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
/*!\brief Checking the capacity of the given custom matrix.
//
// \param matrix The custom matrix to be checked.
// \param minCapacity The expected minimum capacity of the custom matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given custom matrix. In case the actual capacity
// is smaller than the given expected minimum capacity, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >  // Type of the custom matrix
void AlignedUnpaddedTest::checkCapacity( const Type& matrix, size_t minCapacity ) const
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
/*!\brief Checking the number of non-zero elements of the given custom matrix.
//
// \param matrix The custom matrix to be checked.
// \param expectedNonZeros The expected number of non-zero elements of the custom matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given custom matrix. In
// case the actual number of non-zero elements does not correspond to the given expected
// number, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the custom matrix
void AlignedUnpaddedTest::checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const
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
/*!\brief Checking the number of non-zero elements in a specific row/column of the given custom matrix.
//
// \param matrix The custom matrix to be checked.
// \param index The row/column to be checked.
// \param expectedNonZeros The expected number of non-zero elements in the specified row/column.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements in the specified row/column of the
// given custom matrix. In case the actual number of non-zero elements does not correspond
// to the given expected number, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the custom matrix
void AlignedUnpaddedTest::checkNonZeros( const Type& matrix, size_t index, size_t expectedNonZeros ) const
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
/*!\brief Testing the functionality of the aligned/unpadded CustomMatrix class template.
//
// \return void
*/
void runTest()
{
   AlignedUnpaddedTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the aligned/unpadded CustomMatrix class test.
*/
#define RUN_CUSTOMMATRIX_ALIGNED_UNPADDED_TEST \
   blazetest::mathtest::matrices::custommatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace custommatrix

} // namespace matrices

} // namespace mathtest

} // namespace blazetest

#endif
