//=================================================================================================
/*!
//  \file blazetest/mathtest/views/elements/SparseTest.h
//  \brief Header file for the Elements sparse test
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

#ifndef _BLAZETEST_MATHTEST_VIEWS_ELEMENTS_SPARSETEST_H_
#define _BLAZETEST_MATHTEST_VIEWS_ELEMENTS_SPARSETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/Elements.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

namespace views {

namespace elements {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the sparse Elements specialization.
//
// This class represents a test suite for the blaze::Elements class template specialization for
// sparse element selections. It performs a series of both compile time as well as runtime tests.
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
   void testReserve();
   void testSet();
   void testInsert();
   void testAppend();
   void testErase();
   void testFind();
   void testLowerBound();
   void testUpperBound();
   void testIsDefault();
   void testIsSame();
   void testSubvector();
   void testElements();

   template< typename Type >
   void checkSize( const Type& vector, size_t expectedSize ) const;

   template< typename Type >
   void checkCapacity( const Type& vector, size_t minCapacity ) const;

   template< typename Type >
   void checkNonZeros( const Type& vector, size_t expectedNonZeros ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initialize();
   //@}
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   using VT = blaze::CompressedVector<int,blaze::rowVector>;  //!< Compressed row vector type
   using ET = blaze::Elements<VT>;                            //!< Elements type for compressed row vectors.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   VT vec_;  //!< Compressed column vector.
             /*!< The 8-dimensional sparse vector is initialized as
                  \f[\left(\begin{array}{*{4}{c}}
                  0 & 1 & 0 & -2 & -3 & 0 & 4 & 0 \\
                  \end{array}\right)\f]. */

   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( ET );
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
/*!\brief Checking the size of the given sparse vector.
//
// \param vector The sparse vector to be checked.
// \param expectedSize The expected size of the sparse vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the size of the given sparse vector. In case the actual size does not
// correspond to the given expected size, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the sparse vector
void SparseTest::checkSize( const Type& vector, size_t expectedSize ) const
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
/*!\brief Checking the capacity of the given sparse vector.
//
// \param vector The sparse vector to be checked.
// \param minCapacity The expected minimum capacity.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given sparse vector. In case the actual capacity is
// smaller than the given expected minimum capacity, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the sparse vector
void SparseTest::checkCapacity( const Type& vector, size_t minCapacity ) const
{
   if( capacity( vector ) < minCapacity ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Capacity                 : " << capacity( vector ) << "\n"
          << "   Expected minimum capacity: " << minCapacity << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of non-zero elements of the given sparse vector.
//
// \param object The sparse vector to be checked.
// \param expectedNonZeros The expected number of non-zero elements.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given sparse vector. In case
// the actual number of non-zero elements does not correspond to the given expected number,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the sparse vector
void SparseTest::checkNonZeros( const Type& vector, size_t expectedNonZeros ) const
{
   if( nonZeros( vector ) != expectedNonZeros ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of non-zero elements\n"
          << " Details:\n"
          << "   Number of non-zeros         : " << nonZeros( vector ) << "\n"
          << "   Expected number of non-zeros: " << expectedNonZeros << "\n";
      throw std::runtime_error( oss.str() );
   }

   if( capacity( vector ) < nonZeros( vector ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Number of non-zeros: " << nonZeros( vector ) << "\n"
          << "   Capacity           : " << capacity( vector ) << "\n";
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
/*!\brief Testing the functionality of the sparse Elements specialization.
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
/*!\brief Macro for the execution of the Elements sparse test.
*/
#define RUN_ELEMENTS_SPARSE_TEST \
   blazetest::mathtest::views::elements::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace elements

} // namespace views

} // namespace mathtest

} // namespace blazetest

#endif
