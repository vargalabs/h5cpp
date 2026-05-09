//=================================================================================================
/*!
//  \file blazetest/mathtest/matrices/compressedmatrix/ProxyTest.h
//  \brief Header file for the CompressedMatrix proxy test
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

#ifndef _BLAZETEST_MATHTEST_MATRICES_COMPRESSEDMATRIX_PROXYTEST_H_
#define _BLAZETEST_MATHTEST_MATRICES_COMPRESSEDMATRIX_PROXYTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace compressedmatrix {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the access proxy of the CompressedMatrix class template.
//
// This class represents a test suite for the access proxy of the blaze::CompressedMatrix class
// template, the blaze::MatrixAccessProxy class template. It performs a series of both compile
// time as well as runtime tests.
*/
class ProxyTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit ProxyTest();
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
   void testAssignment();
   void testAddAssign();
   void testSubAssign();
   void testMultAssign();
   void testDivAssign();
   void testModAssign();
   void testScaling();
   void testSubscript();
   void testFunctionCall();
   void testIterator();
   void testNonZeros();
   void testReset();
   void testClear();
   void testResize();
   void testExtend();
   void testReserve();
   void testTrim();
   void testSwap();
   void testSet();
   void testInsert();
   void testAppend();
   void testErase();
   void testFind();
   void testLowerBound();
   void testUpperBound();
   void testTranspose();
   void testCTranspose();
   void testInvert();

   template< typename Type >
   void checkSize( const Type& vector, size_t expectedSize ) const;

   template< typename Type >
   void checkRows( const Type& matrix, size_t expectedRows ) const;

   template< typename Type >
   void checkColumns( const Type& matrix, size_t expectedColumns ) const;

   template< typename Type >
   void checkCapacity( const Type& object, size_t minCapacity ) const;

   template< typename Type >
   void checkCapacity( const Type& matrix, size_t index, size_t minCapacity ) const;

   template< typename Type >
   void checkNonZeros( const Type& object, size_t expectedNonZeros ) const;

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
   using DV = blaze::DynamicVector<int,blaze::rowVector>;     //!< Type of the dense vector elements.
   using SV = blaze::CompressedVector<int,blaze::rowVector>;  //!< Type of the sparse vector elements.
   using DM = blaze::DynamicMatrix<int,blaze::rowMajor>;      //!< Type of the dense matrix elements.
   using SM = blaze::CompressedMatrix<int,blaze::rowMajor>;   //!< Type of the sparse matrix elements.

   //! Row-major compressed matrix type with dense vector elements.
   using DVM = blaze::CompressedMatrix<DV,blaze::rowMajor>;

   //! Column-major compressed matrix type with dense vector elements.
   using ODVM = DVM::OppositeType;

   //! Row-major compressed matrix type with sparse vector elements.
   using SVM = blaze::CompressedMatrix<SV,blaze::rowMajor>;

   //! Column-major compressed matrix type with sparse vector elements.
   using OSVM = SVM::OppositeType;

   //! Row-major compressed matrix type with dense matrix elements.
   using DMM = blaze::CompressedMatrix<DM,blaze::rowMajor>;

   //! Column-major compressed matrix type with dense matrix elements.
   using ODMM = DMM::OppositeType;

   //! Row-major compressed matrix type with sparse matrix elements.
   using SMM = blaze::CompressedMatrix<SM,blaze::rowMajor>;

   //! Column-major compressed matrix type with sparse matrix elements.
   using OSMM = SMM::OppositeType;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( DVM  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ODVM );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( SVM  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OSVM );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( DMM  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ODMM );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( SMM  );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OSMM );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( DVM );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( SVM );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( DMM );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( SMM );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ODVM );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OSVM );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ODMM );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OSMM );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DVM, ODVM::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SVM, OSVM::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DMM, ODMM::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SMM, OSMM::TransposeType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DVM::ElementType, ODVM::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SVM::ElementType, OSVM::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( DMM::ElementType, ODMM::ElementType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( SMM::ElementType, OSMM::ElementType );
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
/*!\brief Checking the size of the given vector.
//
// \param vector The vector to be checked.
// \param expectedSize The expected size of the vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the size of the given vector. In case the actual size does not correspond
// to the given expected size, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the vector
void ProxyTest::checkSize( const Type& vector, size_t expectedSize ) const
{
   if( size( vector ) != expectedSize ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid size detected\n"
          << " Details:\n"
          << "   Size         : " << size( vector ) << "\n"
          << "   Expected size: " << expectedSize << "\n";
      throw std::runtime_error( oss.str() );
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
void ProxyTest::checkRows( const Type& matrix, size_t expectedRows ) const
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
template< typename Type >  // Type of the  matrix
void ProxyTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
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
/*!\brief Checking the capacity of the given vector/matrix.
//
// \param object The vector/matrix to be checked.
// \param minCapacity The expected minimum capacity of the vector/matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given vector/matrix. In case the actual capacity is
// smaller than the given expected minimum capacity, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the vector/matrix
void ProxyTest::checkCapacity( const Type& object, size_t minCapacity ) const
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
void ProxyTest::checkCapacity( const Type& matrix, size_t index, size_t minCapacity ) const
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
/*!\brief Checking the number of non-zero elements of the given vector/matrix.
//
// \param object The vector/matrix to be checked.
// \param expectedNonZeros The expected number of non-zero elements of the vector/matrix.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given vector/matrix. In case
// the actual number of non-zero elements does not correspond to the given expected number,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the vector/matrix
void ProxyTest::checkNonZeros( const Type& object, size_t expectedNonZeros ) const
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
/*!\brief Checking the number of non-zero elements in a specific row/column of the given matrix.
//
// \param matrix The matrix to be checked.
// \param index The row/column to be checked.
// \param expectedNonZeros The expected number of non-zero elements in the specified row/column.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements in the specified row/column of the given
// matrix. In case the actual number of non-zero elements does not correspond to the given expected
// number, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the matrix
void ProxyTest::checkNonZeros( const Type& matrix, size_t index, size_t expectedNonZeros ) const
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
/*!\brief Testing the functionality of the MatrixAccessProxy class template.
//
// \return void
*/
void runTest()
{
   ProxyTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the MatrixAccessProxy class test.
*/
#define RUN_COMPRESSEDMATRIX_PROXY_TEST \
   blazetest::mathtest::matrices::compressedmatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace compressedmatrix

} // namespace matrices

} // namespace mathtest

} // namespace blazetest

#endif
