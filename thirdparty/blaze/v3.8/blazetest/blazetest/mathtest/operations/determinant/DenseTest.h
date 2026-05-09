//=================================================================================================
/*!
//  \file blazetest/mathtest/operations/determinant/DenseTest.h
//  \brief Header file for the dense matrix determinant test
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

#ifndef _BLAZETEST_MATHTEST_OPERATIONS_DETERMINANT_DENSETEST_H_
#define _BLAZETEST_MATHTEST_OPERATIONS_DETERMINANT_DENSETEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/Aliases.h>
#include <blaze/math/shims/Equal.h>
#include <blazetest/system/LAPACK.h>


namespace blazetest {

namespace mathtest {

namespace operations {

namespace determinant {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for all dense matrix determinant tests.
//
// This class represents a test suite for the dense matrix determinant functionality. It computes
// the determinant for different matrices using different approaches.
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
   void testSpecific();

   template< typename Type >
   void testRandom2x2();

   template< typename Type >
   void testRandom3x3();

   template< typename Type >
   void testRandom4x4();

   template< typename Type >
   void testRandom5x5();

   template< typename Type >
   void testRandom6x6();
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
/*!\brief Test of the determinant functionality with random 2x2 matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function computes determinants for random 2x2 matrices. In case an error is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom2x2()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::equal;

   using ET = blaze::ElementType_t<Type>;

   Type A;
   resize( A, 2UL, 2UL );
   randomize( A );

   const ET res1( det   ( A ) );
   const ET res2( det2x2( A ) );
   const ET res3( detNxN( A ) );

   if( !equal( res1, res2 ) || !equal( res1, res3 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid determinant evaluation\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( ET ).name() << "\n"
          << "   Result det   (): " << res1 << "\n"
          << "   Result det2x2(): " << res2 << "\n"
          << "   Result detNxN(): " << res3 << "\n";
      throw std::runtime_error( oss.str() );
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the determinant functionality with random 3x3 matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function computes determinants for random 3x3 matrices. In case an error is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom3x3()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::equal;

   using ET = blaze::ElementType_t<Type>;

   Type A;
   resize( A, 3UL, 3UL );
   randomize( A );

   const ET res1( det   ( A ) );
   const ET res2( det3x3( A ) );
   const ET res3( detNxN( A ) );

   if( !equal( res1, res2 ) || !equal( res1, res3 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid determinant evaluation\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( ET ).name() << "\n"
          << "   Result det   (): " << res1 << "\n"
          << "   Result det3x3(): " << res2 << "\n"
          << "   Result detNxN(): " << res3 << "\n";
      throw std::runtime_error( oss.str() );
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the determinant functionality with random 4x4 matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function computes determinants for random 4x4 matrices. In case an error is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom4x4()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::equal;

   using ET = blaze::ElementType_t<Type>;

   Type A;
   resize( A, 4UL, 4UL );
   randomize( A );

   const ET res1( det   ( A ) );
   const ET res2( det4x4( A ) );
   const ET res3( detNxN( A ) );

   if( !equal( res1, res2 ) || !equal( res1, res3 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid determinant evaluation\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( ET ).name() << "\n"
          << "   Result det   (): " << res1 << "\n"
          << "   Result det4x4(): " << res2 << "\n"
          << "   Result detNxN(): " << res3 << "\n";
      throw std::runtime_error( oss.str() );
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the determinant functionality with random 5x5 matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function computes determinants for random 5x5 matrices. In case an error is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom5x5()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::equal;

   using ET = blaze::ElementType_t<Type>;

   Type A;
   resize( A, 5UL, 5UL );
   randomize( A );

   const ET res1( det   ( A ) );
   const ET res2( det5x5( A ) );
   const ET res3( detNxN( A ) );

   if( !equal( res1, res2 ) || !equal( res1, res3 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid determinant evaluation\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( ET ).name() << "\n"
          << "   Result det   (): " << res1 << "\n"
          << "   Result det5x5(): " << res2 << "\n"
          << "   Result detNxN(): " << res3 << "\n";
      throw std::runtime_error( oss.str() );
   }

#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the determinant functionality with random 6x6 matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function computes determinants for random 6x6 matrices. In case an error is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename Type >
void DenseTest::testRandom6x6()
{
#if BLAZETEST_MATHTEST_LAPACK_MODE

   using blaze::equal;

   using ET = blaze::ElementType_t<Type>;

   Type A;
   resize( A, 6UL, 6UL );
   randomize( A );

   const ET res1( det   ( A ) );
   const ET res2( det6x6( A ) );
   const ET res3( detNxN( A ) );

   if( !equal( res1, res2 ) || !equal( res1, res3 ) ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Invalid determinant evaluation\n"
          << " Details:\n"
          << "   Matrix type:\n"
          << "     " << typeid( Type ).name() << "\n"
          << "   Element type:\n"
          << "     " << typeid( ET ).name() << "\n"
          << "   Result det   (): " << res1 << "\n"
          << "   Result det6x6(): " << res2 << "\n"
          << "   Result detNxN(): " << res3 << "\n";
      throw std::runtime_error( oss.str() );
   }

#endif
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the dense matrix determinant.
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
/*!\brief Macro for the execution of the dense matrix determinant test.
*/
#define RUN_DETERMINANT_DENSE_TEST \
   blazetest::mathtest::operations::determinant::runTest()
/*! \endcond */
//*************************************************************************************************

} // namespace determinant

} // namespace operations

} // namespace mathtest

} // namespace blazetest

#endif
