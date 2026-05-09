//=================================================================================================
/*!
//  \file blazetest/mathtest/adaptors/symmetricmatrix/ColumnNonScalarTest.h
//  \brief Header file for the SymmetricMatrix column non-scalar test
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

#ifndef _BLAZETEST_MATHTEST_ADAPTORS_SYMMETRICMATRIX_COLUMNNONSCALARTEST_H_
#define _BLAZETEST_MATHTEST_ADAPTORS_SYMMETRICMATRIX_COLUMNNONSCALARTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/Column.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace adaptors {

namespace symmetricmatrix {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for assignment tests to a single column of a non-scalar SymmetricMatrix.
//
// This class performs assignment tests to a single column of a SymmetricMatrix with non-scalar
// element type. It performs a series of both compile time as well as runtime tests.
*/
class ColumnNonScalarTest
{
 private:
   //**Type definitions****************************************************************************
   //! Type of a resizable, non-scalar element.
   using VT = blaze::DynamicVector<int,blaze::rowVector>;

   //! Type of the dense non-scalar symmetric matrix.
   using DST = blaze::SymmetricMatrix< blaze::DynamicMatrix<VT,blaze::rowMajor> >;

   //! Opposite dense non-scalar symmetric matrix type.
   using DOST = DST::OppositeType;

   //! Type of the sparse non-scalar symmetric matrix.
   using SST = blaze::SymmetricMatrix< blaze::CompressedMatrix<VT,blaze::rowMajor> >;

   //! Opposite sparse non-scalar symmetric matrix type.
   using SOST = SST::OppositeType;
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit ColumnNonScalarTest();
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
   template< typename ST > void testAssignment();

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
   template< typename ST > void init( ST& sym );
   inline VT vec( int value );
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
/*!\brief Test of the assignment to columns of a SymmetricMatrix.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the assignment to a single column of a SymmetricMatrix. In case
// an error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename ST >  // Type of the symmetric matrix
void ColumnNonScalarTest::testAssignment()
{
   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   // ( (  1 ) ( -4 ) (  7 ) )      ( (  1 ) (  2 ) (  7 ) )
   // ( ( -4 ) (  2 ) (    ) )  =>  ( (  2 ) (  8 ) (  4 ) )
   // ( (  7 ) (    ) (  3 ) )      ( (  7 ) (  4 ) (  3 ) )
   {
      test_ = "Dense vector assignment test";

      blaze::DynamicVector<VT,blaze::columnVector> tmp( 3UL );
      tmp[0] = vec( 2 );
      tmp[1] = vec( 8 );
      tmp[2] = vec( 4 );

      ST sym;
      init( sym );

      auto col1 = column( sym, 1UL );
      col1 = tmp;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 9UL );

      if( col1[0] != vec( 2 ) || col1[1] != vec( 8 ) || col1[2] != vec( 4 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( ( 2 ) ( 8 ) ( 4 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) || sym(0,2) != vec( 7 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 8 ) || sym(1,2) != vec( 4 ) ||
          sym(2,0) != vec( 7 ) || sym(2,1) != vec( 4 ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) ( 7 ) )\n"
                                     "( ( 2 ) ( 8 ) ( 4 ) )\n"
                                     "( ( 7 ) ( 4 ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   // ( (  1 ) ( -4 ) (  7 ) )      ( (  1 ) (  2 ) (  7 ) )
   // ( ( -4 ) (  2 ) (    ) )  =>  ( (  2 ) (  8 ) (  4 ) )
   // ( (  7 ) (    ) (  3 ) )      ( (  7 ) (  4 ) (  3 ) )
   {
      test_ = "Dense vector assignment test";

      blaze::CompressedVector<VT,blaze::columnVector> tmp( 3UL, 3UL );
      tmp[0] = vec( 2 );
      tmp[1] = vec( 8 );
      tmp[2] = vec( 4 );

      ST sym;
      init( sym );

      auto col1 = column( sym, 1UL );
      col1 = tmp;

      checkRows    ( sym, 3UL );
      checkColumns ( sym, 3UL );
      checkNonZeros( sym, 9UL );

      if( col1[0] != vec( 2 ) || col1[1] != vec( 8 ) || col1[2] != vec( 4 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << col1 << "\n"
             << "   Expected result:\n( ( 2 ) ( 8 ) ( 4 ) )\n";
         throw std::runtime_error( oss.str() );
      }

      if( sym(0,0) != vec( 1 ) || sym(0,1) != vec( 2 ) || sym(0,2) != vec( 7 ) ||
          sym(1,0) != vec( 2 ) || sym(1,1) != vec( 8 ) || sym(1,2) != vec( 4 ) ||
          sym(2,0) != vec( 7 ) || sym(2,1) != vec( 4 ) || sym(2,2) != vec( 3 ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Row access failed\n"
             << " Details:\n"
             << "   Result:\n" << sym << "\n"
             << "   Expected result:\n( ( 1 ) ( 2 ) ( 7 ) )\n"
                                     "( ( 2 ) ( 8 ) ( 4 ) )\n"
                                     "( ( 7 ) ( 4 ) ( 3 ) )\n";
         throw std::runtime_error( oss.str() );
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
void ColumnNonScalarTest::checkRows( const Type& matrix, size_t expectedRows ) const
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
void ColumnNonScalarTest::checkColumns( const Type& matrix, size_t expectedColumns ) const
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
void ColumnNonScalarTest::checkNonZeros( const Type& matrix, size_t expectedNonZeros ) const
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
/*!\brief Initializing the given symmetric matrix.
//
// \return void
//
// This function is called before each test case to initialize the given symmetric matrix.
*/
template< typename ST >
void ColumnNonScalarTest::init( ST& sym )
{
   sym.resize( 3UL );
   sym(0,0) = vec(  1 );
   sym(0,1) = vec( -4 );
   sym(0,2) = vec(  7 );
   sym(1,1) = vec(  2 );
   sym(2,2) = vec(  3 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setup of a vector.
//
// \param value The value of the vector.
// \return The created vector.
// \exception std::runtime_error Error detected.
//
// This function creates a single vector of size 1. The element of the vector is initialized with
// the given integer value.
*/
inline ColumnNonScalarTest::VT ColumnNonScalarTest::vec( int value )
{
   return VT( 1UL, value );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the assignment to a single column of a non-scalar SymmetricMatrix.
//
// \return void
*/
void runTest()
{
   ColumnNonScalarTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the SymmetricMatrix column non-scalar test.
*/
#define RUN_SYMMETRICMATRIX_COLUMNNONSCALAR_TEST \
   blazetest::mathtest::adaptors::symmetricmatrix::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace symmetricmatrix

} // namespace adaptors

} // namespace mathtest

} // namespace blazetest

#endif
