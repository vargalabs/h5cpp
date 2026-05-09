//=================================================================================================
/*!
//  \file blazetest/mathtest/vectors/sparsevector/ZeroTest.h
//  \brief Header file for the zero SparseVector operation test
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

#ifndef _BLAZETEST_MATHTEST_VECTORS_SPARSEVECTOR_ZEROTEST_H_
#define _BLAZETEST_MATHTEST_VECTORS_SPARSEVECTOR_ZEROTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>


namespace blazetest {

namespace mathtest {

namespace vectors {

namespace sparsevector {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all tests of the SparseVector functionality.
//
// This class represents a test suite for the SparseVector functionality contained in the
// <em><blaze/math/sparse/SparseVector.h></em> header file. It performs a series of runtime
// tests with zero vectors.
*/
class ZeroTest
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit ZeroTest();
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
   void testIsNan();
   void testIsInf();
   void testIsFinite();
   void testIsUniform();
   void testIsZero();
   void testMean();
   void testVar();
   void testStdDev();

   template< typename Type >
   void checkSize( const Type& vector, size_t expectedSize ) const;

   template< typename Type >
   void checkCapacity( const Type& vector, size_t minCapacity ) const;

   template< typename Type >
   void checkNonZeros( const Type& vector, size_t nonzeros ) const;
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
/*!\brief Checking the size of the given sparse vector.
//
// \param vector The sparse vector to be checked.
// \param expectedSize The expected size of the sparse vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the size of the given sparse vector. In case the actual size
// does not correspond to the given expected size, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >  // Type of the sparse vector
void ZeroTest::checkSize( const Type& vector, size_t expectedSize ) const
{
   if( vector.size() != expectedSize ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid size detected\n"
          << " Details:\n"
          << "   Size         : " << vector.size() << "\n"
          << "   Expected size: " << expectedSize << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the capacity of the given sparse vector.
//
// \param vector The sparse vector to be checked.
// \param minCapacity The expected minimum capacity of the sparse vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the capacity of the given sparse vector. In case the actual capacity
// is smaller than the given expected minimum capacity, a \a std::runtime_error exception is
// thrown.
*/
template< typename Type >  // Type of the sparse vector
void ZeroTest::checkCapacity( const Type& vector, size_t minCapacity ) const
{
   if( vector.capacity() < minCapacity ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid capacity detected\n"
          << " Details:\n"
          << "   Capacity                 : " << vector.capacity() << "\n"
          << "   Expected minimum capacity: " << minCapacity << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checking the number of non-zero elements of the given sparse vector.
//
// \param vector The sparse vector to be checked.
// \param expectedNonZeros The expected number of non-zero elements of the sparse vector.
// \return void
// \exception std::runtime_error Error detected.
//
// This function checks the number of non-zero elements of the given sparse vector. In
// case the actual number of non-zero elements does not correspond to the given expected
// number, a \a std::runtime_error exception is thrown.
*/
template< typename Type >  // Type of the sparse vector
void ZeroTest::checkNonZeros( const Type& vector, size_t expectedNonZeros ) const
{
   if( vector.nonZeros() != expectedNonZeros ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid number of non-zero elements\n"
          << " Details:\n"
          << "   Number of non-zeros         : " << vector.nonZeros() << "\n"
          << "   Expected number of non-zeros: " << expectedNonZeros << "\n";
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
/*!\brief Testing the functionality of the SparseVector class template.
//
// \return void
*/
void runTest()
{
   ZeroTest();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of the zero SparseVector operation test.
*/
#define RUN_SPARSEVECTOR_ZERO_TEST \
   blazetest::mathtest::vectors::sparsevector::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace sparsevector

} // namespace vectors

} // namespace mathtest

} // namespace blazetest

#endif
