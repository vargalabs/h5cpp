//=================================================================================================
/*!
//  \file blazetest/mathtest/simd/OperationTest.h
//  \brief Header file for the SIMD operation test
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

#ifndef _BLAZETEST_MATHTEST_SIMD_OPERATIONTEST_H_
#define _BLAZETEST_MATHTEST_SIMD_OPERATIONTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <blaze/math/Accuracy.h>
#include <blaze/math/shims/Abs.h>
#include <blaze/math/shims/Acos.h>
#include <blaze/math/shims/Acosh.h>
#include <blaze/math/shims/Asin.h>
#include <blaze/math/shims/Asinh.h>
#include <blaze/math/shims/Atan.h>
#include <blaze/math/shims/Atan2.h>
#include <blaze/math/shims/Atanh.h>
#include <blaze/math/shims/Cbrt.h>
#include <blaze/math/shims/Ceil.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/Cos.h>
#include <blaze/math/shims/Cosh.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/Erf.h>
#include <blaze/math/shims/Erfc.h>
#include <blaze/math/shims/Exp.h>
#include <blaze/math/shims/Exp2.h>
#include <blaze/math/shims/Floor.h>
#include <blaze/math/shims/Hypot.h>
#include <blaze/math/shims/InvCbrt.h>
#include <blaze/math/shims/InvSqrt.h>
#include <blaze/math/shims/Log.h>
#include <blaze/math/shims/Log2.h>
#include <blaze/math/shims/Log10.h>
#include <blaze/math/shims/Pow.h>
#include <blaze/math/shims/Pow2.h>
#include <blaze/math/shims/Pow3.h>
#include <blaze/math/shims/Pow4.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Round.h>
#include <blaze/math/shims/Sign.h>
#include <blaze/math/shims/Sin.h>
#include <blaze/math/shims/Sinh.h>
#include <blaze/math/shims/Sqrt.h>
#include <blaze/math/shims/Tan.h>
#include <blaze/math/shims/Tanh.h>
#include <blaze/math/shims/Trunc.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/typetraits/HasSIMDAbs.h>
#include <blaze/math/typetraits/HasSIMDAcos.h>
#include <blaze/math/typetraits/HasSIMDAcosh.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDAsin.h>
#include <blaze/math/typetraits/HasSIMDAsinh.h>
#include <blaze/math/typetraits/HasSIMDAtan.h>
#include <blaze/math/typetraits/HasSIMDAtan2.h>
#include <blaze/math/typetraits/HasSIMDAtanh.h>
#include <blaze/math/typetraits/HasSIMDBitand.h>
#include <blaze/math/typetraits/HasSIMDBitor.h>
#include <blaze/math/typetraits/HasSIMDBitxor.h>
#include <blaze/math/typetraits/HasSIMDCbrt.h>
#include <blaze/math/typetraits/HasSIMDCeil.h>
#include <blaze/math/typetraits/HasSIMDConj.h>
#include <blaze/math/typetraits/HasSIMDCos.h>
#include <blaze/math/typetraits/HasSIMDCosh.h>
#include <blaze/math/typetraits/HasSIMDDiv.h>
#include <blaze/math/typetraits/HasSIMDEqual.h>
#include <blaze/math/typetraits/HasSIMDErf.h>
#include <blaze/math/typetraits/HasSIMDErfc.h>
#include <blaze/math/typetraits/HasSIMDExp.h>
#include <blaze/math/typetraits/HasSIMDExp2.h>
#include <blaze/math/typetraits/HasSIMDExp10.h>
#include <blaze/math/typetraits/HasSIMDFloor.h>
#include <blaze/math/typetraits/HasSIMDHypot.h>
#include <blaze/math/typetraits/HasSIMDInvCbrt.h>
#include <blaze/math/typetraits/HasSIMDInvSqrt.h>
#include <blaze/math/typetraits/HasSIMDLog.h>
#include <blaze/math/typetraits/HasSIMDLog2.h>
#include <blaze/math/typetraits/HasSIMDLog10.h>
#include <blaze/math/typetraits/HasSIMDMax.h>
#include <blaze/math/typetraits/HasSIMDMin.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDPow.h>
#include <blaze/math/typetraits/HasSIMDRound.h>
#include <blaze/math/typetraits/HasSIMDShiftLI.h>
#include <blaze/math/typetraits/HasSIMDShiftLV.h>
#include <blaze/math/typetraits/HasSIMDShiftRI.h>
#include <blaze/math/typetraits/HasSIMDShiftRV.h>
#include <blaze/math/typetraits/HasSIMDSign.h>
#include <blaze/math/typetraits/HasSIMDSin.h>
#include <blaze/math/typetraits/HasSIMDSinh.h>
#include <blaze/math/typetraits/HasSIMDSqrt.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/HasSIMDTan.h>
#include <blaze/math/typetraits/HasSIMDTanh.h>
#include <blaze/math/typetraits/HasSIMDTrunc.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/Memory.h>
#include <blaze/util/NonCopyable.h>
#include <blaze/util/policies/Deallocate.h>
#include <blaze/util/Random.h>


namespace blazetest {

namespace mathtest {

namespace simd {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class template for the SIMD operation test.
//
// This class template represents the tests of all available SIMD operations for the given
// numeric data type \a T. In these tests both aligned and unaligned load/store operations
// are used.
*/
template< typename T >  // Data type of the SIMD test
class OperationTest : private blaze::NonCopyable
{
 private:
   //**Type definitions****************************************************************************
   using SIMD     = blaze::SIMDTrait<T>;  //!< SIMD trait for the given numeric type.
   using SIMDType = typename SIMD::Type;  //!< SIMD type for the given numeric type.
   //**********************************************************************************************

   //**********************************************************************************************
   static constexpr size_t SIMDSIZE = SIMD::size;  //!< Number of elements in a single SIMD vector.
   static constexpr size_t N = 256UL;              //!< Number of numeric values to be worked on.
   static constexpr size_t NN = N + SIMDSIZE;      //!< Total number of numeric values in each array.
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit OperationTest();
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~OperationTest();
   //@}
   //**********************************************************************************************

 private:
   //**Test functions******************************************************************************
   /*!\name Test functions */
   //@{
   void testStorea        ();
   void testStoreu        ( size_t offset );
   void testStream        ();
   void testSet           ();

   void testEquality      ( blaze::TrueType , blaze::TrueType  );
   void testEquality      ( blaze::TrueType , blaze::FalseType );
   void testEquality      ( blaze::FalseType, blaze::TrueType  );
   void testEquality      ( blaze::FalseType, blaze::FalseType );
   void testInequality    ( blaze::TrueType , blaze::TrueType  );
   void testInequality    ( blaze::TrueType , blaze::FalseType );
   void testInequality    ( blaze::FalseType, blaze::TrueType  );
   void testInequality    ( blaze::FalseType, blaze::FalseType );

   void testAddition      ( blaze::TrueType  );
   void testAddition      ( blaze::FalseType );
   void testSubtraction   ( blaze::TrueType  );
   void testSubtraction   ( blaze::FalseType );
   void testMultiplication( blaze::TrueType  );
   void testMultiplication( blaze::FalseType );
   void testFmadd         ( blaze::TrueType  );
   void testFmadd         ( blaze::FalseType );
   void testFmsub         ( blaze::TrueType  );
   void testFmsub         ( blaze::FalseType );
   void testDivision      ( blaze::TrueType  );
   void testDivision      ( blaze::FalseType );

   void testBitand        ( blaze::TrueType  );
   void testBitand        ( blaze::FalseType );
   void testBitor         ( blaze::TrueType  );
   void testBitor         ( blaze::FalseType );
   void testBitxor        ( blaze::TrueType  );
   void testBitxor        ( blaze::FalseType );

   void testShiftLI       ( blaze::TrueType  );
   void testShiftLI       ( blaze::FalseType );
   void testShiftLV       ( blaze::TrueType  );
   void testShiftLV       ( blaze::FalseType );
   void testShiftRI       ( blaze::TrueType  );
   void testShiftRI       ( blaze::FalseType );
   void testShiftRV       ( blaze::TrueType  );
   void testShiftRV       ( blaze::FalseType );

   void testMin           ( blaze::TrueType  );
   void testMin           ( blaze::FalseType );
   void testMax           ( blaze::TrueType  );
   void testMax           ( blaze::FalseType );

   void testAbs           ( blaze::TrueType  );
   void testAbs           ( blaze::FalseType );
   void testSign          ( blaze::TrueType  );
   void testSign          ( blaze::FalseType );

   void testFloor         ( blaze::TrueType  );
   void testFloor         ( blaze::FalseType );
   void testCeil          ( blaze::TrueType  );
   void testCeil          ( blaze::FalseType );
   void testTrunc         ( blaze::TrueType  );
   void testTrunc         ( blaze::FalseType );
   void testRound         ( blaze::TrueType  );
   void testRound         ( blaze::FalseType );

   void testConj          ( blaze::TrueType  );
   void testConj          ( blaze::FalseType );
   void testSqrt          ( blaze::TrueType  );
   void testSqrt          ( blaze::FalseType );
   void testInvSqrt       ( blaze::TrueType  );
   void testInvSqrt       ( blaze::FalseType );
   void testCbrt          ( blaze::TrueType  );
   void testCbrt          ( blaze::FalseType );
   void testInvCbrt       ( blaze::TrueType  );
   void testInvCbrt       ( blaze::FalseType );
   void testHypot         ( blaze::TrueType  );
   void testHypot         ( blaze::FalseType );

   void testPow           ( blaze::TrueType  );
   void testPow           ( blaze::FalseType );
   void testPow2          ( blaze::TrueType  );
   void testPow2          ( blaze::FalseType );
   void testPow3          ( blaze::TrueType  );
   void testPow3          ( blaze::FalseType );
   void testPow4          ( blaze::TrueType  );
   void testPow4          ( blaze::FalseType );

   void testExp           ( blaze::TrueType  );
   void testExp           ( blaze::FalseType );
   void testExp2          ( blaze::TrueType  );
   void testExp2          ( blaze::FalseType );
   void testExp10         ( blaze::TrueType  );
   void testExp10         ( blaze::FalseType );

   void testLog           ( blaze::TrueType  );
   void testLog           ( blaze::FalseType );
   void testLog2          ( blaze::TrueType  );
   void testLog2          ( blaze::FalseType );
   void testLog10         ( blaze::TrueType  );
   void testLog10         ( blaze::FalseType );

   void testSin           ( blaze::TrueType  );
   void testSin           ( blaze::FalseType );
   void testAsin          ( blaze::TrueType  );
   void testAsin          ( blaze::FalseType );
   void testSinh          ( blaze::TrueType  );
   void testSinh          ( blaze::FalseType );
   void testAsinh         ( blaze::TrueType  );
   void testAsinh         ( blaze::FalseType );

   void testCos           ( blaze::TrueType  );
   void testCos           ( blaze::FalseType );
   void testAcos          ( blaze::TrueType  );
   void testAcos          ( blaze::FalseType );
   void testCosh          ( blaze::TrueType  );
   void testCosh          ( blaze::FalseType );
   void testAcosh         ( blaze::TrueType  );
   void testAcosh         ( blaze::FalseType );

   void testTan           ( blaze::TrueType  );
   void testTan           ( blaze::FalseType );
   void testAtan          ( blaze::TrueType  );
   void testAtan          ( blaze::FalseType );
   void testTanh          ( blaze::TrueType  );
   void testTanh          ( blaze::FalseType );
   void testAtanh         ( blaze::TrueType  );
   void testAtanh         ( blaze::FalseType );
   void testAtan2         ( blaze::TrueType  );
   void testAtan2         ( blaze::FalseType );

   void testErf           ( blaze::TrueType  );
   void testErf           ( blaze::FalseType );
   void testErfc          ( blaze::TrueType  );
   void testErfc          ( blaze::FalseType );

   void testSum           ();
   void testProd          ();
   //@}
   //**********************************************************************************************

   //**Error detection functions*******************************************************************
   /*!\name Error detection functions */
   //@{
   void compare( const T* expected, const T* actual ) const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void initialize();
   void initialize( T min, T max );
   void initialize( T* array );
   void initialize( T* array, T min, T max );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   T* a_;  //!< The first aligned array of size NN.
   T* b_;  //!< The second aligned array of size NN.
   T* c_;  //!< The third aligned array of size NN.
   T* d_;  //!< The fourth aligned array of size NN.
   T* e_;  //!< The fifth aligned array of size NN.

   std::string test_;  //!< Label of the currently performed test.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( T );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SIMD operation test.
//
// \exception std::runtime_error Operation error detected.
*/
template< typename T >  // Data type of the SIMD test
OperationTest<T>::OperationTest()
   : a_    ( blaze::allocate<T>( NN ) )  // The first aligned array of size NN
   , b_    ( blaze::allocate<T>( NN ) )  // The second aligned array of size NN
   , c_    ( blaze::allocate<T>( NN ) )  // The third aligned array of size NN
   , d_    ( blaze::allocate<T>( NN ) )  // The fourth aligned array of size NN
   , e_    ( blaze::allocate<T>( NN ) )  // The fifth aligned array of size NN
   , test_ ()                            // Label of the currently performed test
{
   testStorea();

   for( size_t offset=0UL; offset<SIMDSIZE; ++offset ) {
      testStoreu( offset );
   }

   testStream        ();
   testSet           ();

   testEquality      ( blaze::HasSIMDEqual<T,T>(), blaze::IsFloatingPoint<T>() );
   testInequality    ( blaze::HasSIMDEqual<T,T>(), blaze::IsFloatingPoint<T>() );

   testAddition      ( blaze::HasSIMDAdd <T,T>() );
   testSubtraction   ( blaze::HasSIMDSub <T,T>() );
   testMultiplication( blaze::HasSIMDMult<T,T>() );
   testFmadd         ( blaze::HasSIMDMult<T,T>() );
   testFmsub         ( blaze::HasSIMDMult<T,T>() );
   testDivision      ( blaze::HasSIMDDiv <T,T>() );

   testBitand        ( blaze::HasSIMDBitand<T,T>() );
   testBitor         ( blaze::HasSIMDBitor <T,T>() );
   testBitxor        ( blaze::HasSIMDBitxor<T,T>() );

   testShiftLI       ( blaze::HasSIMDShiftLI< T >() );
   testShiftLV       ( blaze::HasSIMDShiftLV<T,T>() );
   testShiftRI       ( blaze::HasSIMDShiftRI< T >() );
   testShiftRV       ( blaze::HasSIMDShiftRV<T,T>() );

   testMin           ( blaze::HasSIMDMin<T,T>() );
   testMax           ( blaze::HasSIMDMax<T,T>() );

   testAbs           ( blaze::HasSIMDAbs < T >() );
   testSign          ( blaze::HasSIMDSign< T >() );

   testFloor         ( blaze::HasSIMDFloor< T >() );
   testCeil          ( blaze::HasSIMDCeil < T >() );
   testTrunc         ( blaze::HasSIMDTrunc< T >() );
   testRound         ( blaze::HasSIMDRound< T >() );

   testConj          ( blaze::HasSIMDConj   < T >() );
   testSqrt          ( blaze::HasSIMDSqrt   < T >() );
   testInvSqrt       ( blaze::HasSIMDInvSqrt< T >() );
   testCbrt          ( blaze::HasSIMDCbrt   < T >() );
   testInvCbrt       ( blaze::HasSIMDInvCbrt< T >() );
   testHypot         ( blaze::HasSIMDHypot  <T,T>() );

   testPow           ( blaze::HasSIMDPow <T,T>() );
   testPow2          ( blaze::HasSIMDMult<T,T>() );
   testPow3          ( blaze::HasSIMDMult<T,T>() );
   testPow4          ( blaze::HasSIMDMult<T,T>() );

   testExp           ( blaze::HasSIMDExp  < T >() );
   testExp2          ( blaze::HasSIMDExp2 < T >() );
   testExp10         ( blaze::HasSIMDExp10< T >() );
   testLog           ( blaze::HasSIMDLog  < T >() );
   testLog2          ( blaze::HasSIMDLog2 < T >() );
   testLog10         ( blaze::HasSIMDLog10< T >() );

   testSin           ( blaze::HasSIMDSin  < T >() );
   testAsin          ( blaze::HasSIMDAsin < T >() );
   testSinh          ( blaze::HasSIMDSinh < T >() );
   testAsinh         ( blaze::HasSIMDAsinh< T >() );

   testCos           ( blaze::HasSIMDCos  < T >() );
   testAcos          ( blaze::HasSIMDAcos < T >() );
   testCosh          ( blaze::HasSIMDCosh < T >() );
   testAcosh         ( blaze::HasSIMDAcosh< T >() );

   testTan           ( blaze::HasSIMDTan  < T >() );
   testAtan          ( blaze::HasSIMDAtan < T >() );
   testTanh          ( blaze::HasSIMDTanh < T >() );
   testAtanh         ( blaze::HasSIMDAtanh< T >() );
   testAtan2         ( blaze::HasSIMDAtan2<T,T>() );

   testErf           ( blaze::HasSIMDErf < T >() );
   testErfc          ( blaze::HasSIMDErfc< T >() );

   testSum           ();
   testProd          ();
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for the SIMD operation test.
*/
template< typename T >  // Data type of the SIMD test
OperationTest<T>::~OperationTest()
{
   blaze::deallocate( a_ );
   blaze::deallocate( b_ );
   blaze::deallocate( c_ );
   blaze::deallocate( d_ );
   blaze::deallocate( e_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the aligned store operation.
//
// \return void
// \exception std::runtime_error Load/store error detected.
//
// This function tests the aligned store operation by copying one array to another via aligned
// load and store. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testStorea()
{
   using blaze::loada;
   using blaze::storea;

   test_  = "storea() operation";

   initialize();

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( b_+i, loada( a_+i ) );
   }

   compare( a_, b_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the aligned, non-temporal store operation.
//
// \return void
// \exception std::runtime_error Load/store error detected.
//
// This function tests the aligned, non-temporal store operation by copying one array to another
// via aligned load and non-temporal store. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testStream()
{
   using blaze::loada;
   using blaze::stream;

   test_  = "stream() operation";

   initialize();

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      stream( b_+i, loada( a_+i ) );
   }

   compare( a_, b_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the unaligned store operation.
//
// \return void
// \exception std::runtime_error Load/store error detected.
//
// This function tests the unaligned store operation by copying one array to another via unaligned
// load and store. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testStoreu( size_t offset )
{
   using blaze::loadu;
   using blaze::storeu;

   test_  = "storeu() operation";

   initialize();

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storeu( b_+offset+i, loadu( a_+offset+i ) );
   }

   compare( a_+offset, b_+offset );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the set operation.
//
// \return void
// \exception std::runtime_error Set error detected.
//
// This function tests the set operation by comparing the results of a vectorized and a scalar
// array assignment. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSet()
{
   using blaze::set;
   using blaze::storea;
   using blaze::rand;

   test_  = "set() operation";

   initialize();

   const T value( rand<T>() );

   for( size_t i=0UL; i<N; ++i ) {
      b_[i] = value;
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( c_+i, set( value ) );
   }

   compare( b_, c_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the equality comparison.
//
// \return void
// \exception std::runtime_error Comparison error detected.
//
// This function tests the equality comparison for the given floating point data type \a T. In
// case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testEquality( blaze::TrueType, blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Equality comparison";

   {
      initialize();

      for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
         if( !( loada( a_+i ) == loada( a_+i ) ) ) {
            std::ostringstream oss;
            oss.precision( 20 );
            oss << " Test : " << test_ << "\n"
                << " Error: Equality comparison failed at index " << i << "\n"
                << " Details:\n"
                << "   a[" << i << "] = " << a_[i] << "\n"
                << "   b[" << i << "] = " << a_[i] << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
         if( loada( a_+i ) == loada( b_+i ) ) {
            std::ostringstream oss;
            oss.precision( 20 );
            oss << " Test : " << test_ << "\n"
                << " Error: Equality comparison failed at index " << i << "\n"
                << " Details:\n"
                << "   a[" << i << "] = " << a_[i] << "\n"
                << "   b[" << i << "] = " << b_[i] << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
         if( !equal( loada( a_+i ), loada( a_+i ) ) ) {
            std::ostringstream oss;
            oss.precision( 20 );
            oss << " Test : " << test_ << "\n"
                << " Error: Equality comparison failed at index " << i << "\n"
                << " Details:\n"
                << "   a[" << i << "] = " << a_[i] << "\n"
                << "   b[" << i << "] = " << a_[i] << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
         if( equal( loada( a_+i ), loada( b_+i ) ) ) {
            std::ostringstream oss;
            oss.precision( 20 );
            oss << " Test : " << test_ << "\n"
                << " Error: Equality comparison failed at index " << i << "\n"
                << " Details:\n"
                << "   a[" << i << "] = " << a_[i] << "\n"
                << "   b[" << i << "] = " << b_[i] << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      const T accu( blaze::accuracy );
      const T half( 0.5 );

      initialize( -half*accu, half*accu );

      for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
         if( !( loada( a_+i ) == loada( a_+i ) ) ) {
            std::ostringstream oss;
            oss.precision( 20 );
            oss << " Test : " << test_ << "\n"
                << " Error: Equality comparison failed at index " << i << "\n"
                << " Details:\n"
                << "   a[" << i << "] = " << a_[i] << "\n"
                << "   b[" << i << "] = " << a_[i] << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
         if( loada( a_+i ) == loada( b_+i ) ) {
            std::ostringstream oss;
            oss.precision( 20 );
            oss << " Test : " << test_ << "\n"
                << " Error: Equality comparison failed at index " << i << "\n"
                << " Details:\n"
                << "   a[" << i << "] = " << a_[i] << "\n"
                << "   b[" << i << "] = " << b_[i] << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
         if( !equal( loada( a_+i ), loada( a_+i ) ) ) {
            std::ostringstream oss;
            oss.precision( 20 );
            oss << " Test : " << test_ << "\n"
                << " Error: Equality comparison failed at index " << i << "\n"
                << " Details:\n"
                << "   a[" << i << "] = " << a_[i] << "\n"
                << "   b[" << i << "] = " << a_[i] << "\n";
            throw std::runtime_error( oss.str() );
         }
      }

      for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
         if( !equal( loada( a_+i ), loada( b_+i ) ) ) {
            std::ostringstream oss;
            oss.precision( 20 );
            oss << " Test : " << test_ << "\n"
                << " Error: Equality comparison failed at index " << i << "\n"
                << " Details:\n"
                << "   a[" << i << "] = " << a_[i] << "\n"
                << "   b[" << i << "] = " << b_[i] << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the equality comparison.
//
// \return void
// \exception std::runtime_error Comparison error detected.
//
// This function tests the equality comparison for the given integral data type \a T. In case
// any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testEquality( blaze::TrueType, blaze::FalseType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Equality comparison";

   initialize();

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      if( !( loada( a_+i ) == loada( a_+i ) ) ) {
         std::ostringstream oss;
         oss.precision( 20 );
         oss << " Test : " << test_ << "\n"
             << " Error: Equality comparison failed at index " << i << "\n"
             << " Details:\n"
             << "   a[" << i << "] = " << a_[i] << "\n"
             << "   b[" << i << "] = " << a_[i] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      if( loada( a_+i ) == loada( b_+i ) ) {
         std::ostringstream oss;
         oss.precision( 20 );
         oss << " Test : " << test_ << "\n"
             << " Error: Equality comparison failed at index " << i << "\n"
             << " Details:\n"
             << "   a[" << i << "] = " << a_[i] << "\n"
             << "   b[" << i << "] = " << b_[i] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      if( !( equal( loada( a_+i ), loada( a_+i ) ) ) ) {
         std::ostringstream oss;
         oss.precision( 20 );
         oss << " Test : " << test_ << "\n"
             << " Error: Equality comparison failed at index " << i << "\n"
             << " Details:\n"
             << "   a[" << i << "] = " << a_[i] << "\n"
             << "   b[" << i << "] = " << a_[i] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      if( equal( loada( a_+i ), loada( b_+i ) ) ) {
         std::ostringstream oss;
         oss.precision( 20 );
         oss << " Test : " << test_ << "\n"
             << " Error: Equality comparison failed at index " << i << "\n"
             << " Details:\n"
             << "   a[" << i << "] = " << a_[i] << "\n"
             << "   b[" << i << "] = " << b_[i] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the equality comparison.
//
// \return void
//
// This function is called in case the equality comparison is not available for the given
// floating point data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testEquality( blaze::FalseType, blaze::TrueType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the equality comparison.
//
// \return void
//
// This function is called in case the equality comparison is not available for the given
// integral data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testEquality( blaze::FalseType, blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inequality comparison.
//
// \return void
// \exception std::runtime_error Comparison error detected.
//
// This function tests the inequality comparison for the given floating point data type \a T. In
// case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInequality( blaze::TrueType, blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Inequality comparison";

   initialize();

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      if( loada( a_+i ) != loada( a_+i ) ) {
         std::ostringstream oss;
         oss.precision( 20 );
         oss << " Test : " << test_ << "\n"
             << " Error: Inequality comparison failed at index " << i << "\n"
             << " Details:\n"
             << "   a[" << i << "] = " << a_[i] << "\n"
             << "   b[" << i << "] = " << a_[i] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      if( !( loada( a_+i ) != loada( b_+i ) ) ) {
         std::ostringstream oss;
         oss.precision( 20 );
         oss << " Test : " << test_ << "\n"
             << " Error: Inequality comparison failed at index " << i << "\n"
             << " Details:\n"
             << "   a[" << i << "] = " << a_[i] << "\n"
             << "   b[" << i << "] = " << b_[i] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inequality comparison.
//
// \return void
// \exception std::runtime_error Comparison error detected.
//
// This function tests the inequality comparison for the given integral data type \a T. In case
// any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInequality( blaze::TrueType, blaze::FalseType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Inequality comparison";

   initialize();

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      if( loada( a_+i ) != loada( a_+i ) ) {
         std::ostringstream oss;
         oss.precision( 20 );
         oss << " Test : " << test_ << "\n"
             << " Error: Inequality comparison failed at index " << i << "\n"
             << " Details:\n"
             << "   a[" << i << "] = " << a_[i] << "\n"
             << "   b[" << i << "] = " << a_[i] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      if( !( loada( a_+i ) != loada( b_+i ) ) ) {
         std::ostringstream oss;
         oss.precision( 20 );
         oss << " Test : " << test_ << "\n"
             << " Error: Inequality comparison failed at index " << i << "\n"
             << " Details:\n"
             << "   a[" << i << "] = " << a_[i] << "\n"
             << "   b[" << i << "] = " << b_[i] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inequality comparison.
//
// \return void
//
// This function is called in case the inequality comparison is not available for the given
// floating point data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInequality( blaze::FalseType, blaze::TrueType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inequality comparison.
//
// \return void
//
// This function is called in case the inequality comparison is not available for the given
// integral data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInequality( blaze::FalseType, blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the addition operation.
//
// \return void
// \exception std::runtime_error Addition error detected.
//
// This function tests the addition operation by comparing the results of a vectorized and a
// scalar addition. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAddition( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Addition operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] + b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) + loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the addition operation.
//
// \return void
//
// This function is called in case the addition operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAddition( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the subtraction operation.
//
// \return void
// \exception std::runtime_error Subtraction error detected.
//
// This function tests the subtraction operation by comparing the results of a vectorized and a
// scalar subtraction. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSubtraction( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Subtraction operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] - b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) - loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the subtraction operation.
//
// \return void
//
// This function is called in case the subtraction operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSubtraction( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the multiplication operation.
//
// \return void
// \exception std::runtime_error Multiplication error detected.
//
// This function tests the multiplication operation by comparing the results of a vectorized and
// a scalar multiplication. In case any error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testMultiplication( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Multiplication operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] * b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) * loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the multiplication operation.
//
// \return void
//
// This function is called in case the multiplication operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testMultiplication( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the fused multiply-add operation.
//
// \return void
// \exception std::runtime_error Fused multiply-add error detected.
//
// This function tests the fused multiply-add operation by comparing the results of a vectorized
// and a scalar operation. In case any error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testFmadd( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Fused multiply-add operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      d_[i] = a_[i] * b_[i] + c_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( e_+i, loada( a_+i ) * loada( b_+i ) + loada( c_+i ) );
   }

   compare( d_, e_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the fused multiply-add operation.
//
// \return void
//
// This function is called in case the fused multiply-add operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testFmadd( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the fused multiply-subtract operation.
//
// \return void
// \exception std::runtime_error Fused multiply-subtract error detected.
//
// This function tests the fused multiply-subtract operation by comparing the results of a
// vectorized and a scalar operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testFmsub( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Fused multiply-subtract operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      d_[i] = a_[i] * b_[i] - c_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( e_+i, loada( a_+i ) * loada( b_+i ) - loada( c_+i ) );
   }

   compare( d_, e_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the fused multiply-subtract operation.
//
// \return void
//
// This function is called in case the fused multiply-subtract operation is not available for
// the given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testFmsub( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the division operation.
//
// \return void
// \exception std::runtime_error Division error detected.
//
// This function tests the division operation by comparing the results of a vectorized and a
// scalar division. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testDivision( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Division operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] / b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) / loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the division operation.
//
// \return void
//
// This function is called in case the division operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testDivision( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the bitwise AND ('&') operation.
//
// \return void
// \exception std::runtime_error Error in AND computation detected.
//
// This function tests the bitwise AND ('&') operation by comparing the results of a vectorized
// and a scalar operation. In case any error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testBitand( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Bitwise AND ('&') operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] & b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) & loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the bitwise AND ('&') operation.
//
// \return void
//
// This function is called in case the bitwise AND ('&') operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testBitand( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the bitwise OR ('|') operation.
//
// \return void
// \exception std::runtime_error Error in OR computation detected.
//
// This function tests the bitwise OR ('|') operation by comparing the results of a vectorized
// and a scalar operation. In case any error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testBitor( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Bitwise OR ('|') operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] | b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) | loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the bitwise OR ('|') operation.
//
// \return void
//
// This function is called in case the bitwise OR ('|') operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testBitor( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the bitwise XOR ('^') operation.
//
// \return void
// \exception std::runtime_error Error in XOR computation detected.
//
// This function tests the bitwise XOR ('^') operation by comparing the results of a vectorized
// and a scalar operation. In case any error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testBitxor( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Bitwise XOR ('^') operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] ^ b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) ^ loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the bitwise XOR ('^') operation.
//
// \return void
//
// This function is called in case the bitwise XOR ('^') operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testBitxor( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the uniform left-shift operation.
//
// \return void
// \exception std::runtime_error Error in left-shift operation detected.
//
// This function tests the uniform left-shift operation by comparing the results of a vectorized
// and a scalar operation. In case any error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testShiftLI( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Uniform left-shift operation";

   initialize();

   const int shift = blaze::rand<int>( 0, sizeof(T)*8UL-1UL );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] << shift;
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) << shift );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the uniform left-shift operation.
//
// \return void
//
// This function is called in case the uniform left-shift operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testShiftLI( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the componentwise left-shift operation.
//
// \return void
// \exception std::runtime_error Error in left-shift operation detected.
//
// This function tests the componentwise left-shift operation by comparing the results of a
// vectorized and a scalar operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testShiftLV( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Componentwise left-shift operation";

   initialize( a_ );
   initialize( b_, 0, sizeof(T)*8UL-1UL );
   initialize( c_ );
   initialize( d_ );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] << b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) << loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the componentwise left-shift operation.
//
// \return void
//
// This function is called in case the componentwise left-shift operation is not available for
// the given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testShiftLV( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the uniform right-shift operation.
//
// \return void
// \exception std::runtime_error Error in right-shift operation detected.
//
// This function tests the uniform right-shift operation by comparing the results of a vectorized
// and a scalar operation. In case any error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testShiftRI( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Uniform right-shift operation";

   initialize();

   const int shift = blaze::rand<int>( 0, sizeof(T)*8UL-1UL );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] >> shift;
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) >> shift );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the uniform right-shift operation.
//
// \return void
//
// This function is called in case the uniform right-shift operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testShiftRI( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the componentwise right-shift operation.
//
// \return void
// \exception std::runtime_error Error in right-shift operation detected.
//
// This function tests the componentwise right-shift operation by comparing the results of a
// vectorized and a scalar operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testShiftRV( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;

   test_ = "Componentwise right-shift operation";

   initialize( a_ );
   initialize( b_, 0, sizeof(T)*8UL-1UL );
   initialize( c_ );
   initialize( d_ );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = a_[i] >> b_[i];
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, loada( a_+i ) >> loada( b_+i ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the componentwise right-shift operation.
//
// \return void
//
// This function is called in case the componentwise right-shift operation is not available for
// the given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testShiftRV( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the mininum operation.
//
// \return void
// \exception std::runtime_error Error in minimum computation detected.
//
// This function tests the minimum operation by comparing the results of a vectorized and a
// scalar minimum operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testMin( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;
   using blaze::min;

   test_ = "Minimum operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = min( a_[i], b_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, min( loada( a_+i ), loada( b_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the minimum operation.
//
// \return void
//
// This function is called in case the minimum operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testMin( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the maxinum operation.
//
// \return void
// \exception std::runtime_error Error in maximum computation detected.
//
// This function tests the maximum operation by comparing the results of a vectorized and a
// scalar maximum operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testMax( blaze::TrueType )
{
   using blaze::loada;
   using blaze::storea;
   using blaze::max;

   test_ = "Maximum operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = max( a_[i], b_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, max( loada( a_+i ), loada( b_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the maximum operation.
//
// \return void
//
// This function is called in case the maximum operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testMax( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the absolute value operation.
//
// \return void
// \exception std::runtime_error Error in absolute value computation detected.
//
// This function tests the absolute value operation by comparing the results of a vectorized
// and a scalar absolute value operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAbs( blaze::TrueType )
{
   using blaze::abs;
   using blaze::loada;
   using blaze::storea;

   test_ = "Absolute value operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = abs( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, abs( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the absolute value operation.
//
// \return void
//
// This function is called in case the absolute value operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAbs( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the sign operation.
//
// \return void
// \exception std::runtime_error Error in sign computation detected.
//
// This function tests the sign operation by comparing the results of a vectorized and a scalar
// sign operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSign( blaze::TrueType )
{
   using blaze::sign;
   using blaze::loada;
   using blaze::storea;

   test_ = "Sign operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = sign( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, sign( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the sign operation.
//
// \return void
//
// This function is called in case the sign operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSign( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the floor operation.
//
// \return void
// \exception std::runtime_error Error in floor computation detected.
//
// This function tests the floor operation by comparing the results of a vectorized and a scalar
// floor operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testFloor( blaze::TrueType )
{
   using blaze::floor;
   using blaze::loada;
   using blaze::storea;

   test_ = "Floor operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = floor( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, floor( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the floor operation.
//
// \return void
//
// This function is called in case the floor operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testFloor( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the ceil operation.
//
// \return void
// \exception std::runtime_error Error in ceil computation detected.
//
// This function tests the ceil operation by comparing the results of a vectorized and a scalar
// ceil operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCeil( blaze::TrueType )
{
   using blaze::ceil;
   using blaze::loada;
   using blaze::storea;

   test_ = "Ceil operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = ceil( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, ceil( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the ceil operation.
//
// \return void
//
// This function is called in case the ceil operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCeil( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the trunc operation.
//
// \return void
// \exception std::runtime_error Error in trunc computation detected.
//
// This function tests the trunc operation by comparing the results of a vectorized and a scalar
// trunc operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testTrunc( blaze::TrueType )
{
   using blaze::trunc;
   using blaze::loada;
   using blaze::storea;

   test_ = "Trunc operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = trunc( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, trunc( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the trunc operation.
//
// \return void
//
// This function is called in case the trunc operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testTrunc( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the round operation.
//
// \return void
// \exception std::runtime_error Error in round computation detected.
//
// This function tests the round operation by comparing the results of a vectorized and a scalar
// round operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testRound( blaze::TrueType )
{
   using blaze::round;
   using blaze::loada;
   using blaze::storea;

   test_ = "Round operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = round( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, round( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the round operation.
//
// \return void
//
// This function is called in case the round operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testRound( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the conjugate operation.
//
// \return void
// \exception std::runtime_error Error in conjugate computation detected.
//
// This function tests the conjugate operation by comparing the results of a vectorized and a
// scalar conjugate operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testConj( blaze::TrueType )
{
   using blaze::conj;
   using blaze::loada;
   using blaze::storea;

   test_ = "Conjugate operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = conj( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, conj( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the conjugate operation.
//
// \return void
//
// This function is called in case the conjugate operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testConj( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the square root operation.
//
// \return void
// \exception std::runtime_error Error in square root computation detected.
//
// This function tests the square root operation by comparing the results of a vectorized and
// a scalar square root operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSqrt( blaze::TrueType )
{
   using blaze::sqrt;
   using blaze::loada;
   using blaze::storea;

   test_ = "Square root operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = sqrt( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, sqrt( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the square root operation.
//
// \return void
//
// This function is called in case the square root operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSqrt( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse square root operation.
//
// \return void
// \exception std::runtime_error Error in inverse square root computation detected.
//
// This function tests the inverse square root operation by comparing the results of a vectorized
// and a scalar inverse square root operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvSqrt( blaze::TrueType )
{
   using blaze::invsqrt;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse square root operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = invsqrt( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, invsqrt( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse square root operation.
//
// \return void
//
// This function is called in case the inverse square root operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvSqrt( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the cubic root operation.
//
// \return void
// \exception std::runtime_error Error in cubic root computation detected.
//
// This function tests the cubic root operation by comparing the results of a vectorized and
// a scalar cubic root operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCbrt( blaze::TrueType )
{
   using blaze::cbrt;
   using blaze::loada;
   using blaze::storea;

   test_ = "Cubic root operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = cbrt( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, cbrt( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the cubic root operation.
//
// \return void
//
// This function is called in case the cubic root operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCbrt( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse cubic root operation.
//
// \return void
// \exception std::runtime_error Error in inverse cubic root computation detected.
//
// This function tests the inverse cubic root operation by comparing the results of a vectorized
// and a scalar inverse cubic root operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvCbrt( blaze::TrueType )
{
   using blaze::invcbrt;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse cubic root operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = invcbrt( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, invcbrt( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse cubic root operation.
//
// \return void
//
// This function is called in case the inverse cubic root operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testInvCbrt( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the hypotenous operation.
//
// \return void
// \exception std::runtime_error Error in hypotenous computation detected.
//
// This function tests the hypotenous operation by comparing the results of a vectorized and a
// scalar hypotenous operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testHypot( blaze::TrueType )
{
   using blaze::hypot;
   using blaze::loada;
   using blaze::storea;

   test_ = "Hypot operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = hypot( a_[i], b_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, hypot( loada( a_+i ), loada( b_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the hypotenous operation.
//
// \return void
//
// This function is called in case the hypotenous operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testHypot( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the power operation.
//
// \return void
// \exception std::runtime_error Error in power computation detected.
//
// This function tests the power operation by comparing the results of a vectorized and a
// scalar power operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testPow( blaze::TrueType )
{
   using blaze::pow;
   using blaze::loada;
   using blaze::storea;

   test_ = "Power operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = pow( a_[i], b_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, pow( loada( a_+i ), loada( b_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the power operation.
//
// \return void
//
// This function is called in case the power operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testPow( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the pow2 operation.
//
// \return void
// \exception std::runtime_error Error in pow2 computation detected.
//
// This function tests the pow2 operation by comparing the results of a vectorized and a
// scalar pow2 operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testPow2( blaze::TrueType )
{
   using blaze::pow2;
   using blaze::loada;
   using blaze::storea;

   test_ = "Pow2 operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = pow2( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, pow2( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the pow2 operation.
//
// \return void
//
// This function is called in case the pow2 operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testPow2( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the pow3 operation.
//
// \return void
// \exception std::runtime_error Error in pow3 computation detected.
//
// This function tests the pow3 operation by comparing the results of a vectorized and a
// scalar pow3 operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testPow3( blaze::TrueType )
{
   using blaze::pow3;
   using blaze::loada;
   using blaze::storea;

   test_ = "Pow3 operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = pow3( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, pow3( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the pow3 operation.
//
// \return void
//
// This function is called in case the pow3 operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testPow3( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the pow4 operation.
//
// \return void
// \exception std::runtime_error Error in pow4 computation detected.
//
// This function tests the pow4 operation by comparing the results of a vectorized and a
// scalar pow4 operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testPow4( blaze::TrueType )
{
   using blaze::pow4;
   using blaze::loada;
   using blaze::storea;

   test_ = "Pow4 operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = pow4( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, pow4( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the pow4 operation.
//
// \return void
//
// This function is called in case the pow4 operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testPow4( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the exp() operation.
//
// \return void
// \exception std::runtime_error Error in exp() computation detected.
//
// This function tests the exp() operation by comparing the results of a vectorized and a
// scalar exp() operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testExp( blaze::TrueType )
{
   using blaze::exp;
   using blaze::loada;
   using blaze::storea;

   test_ = "Exponent operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = exp( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, exp( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the exp() operation.
//
// \return void
//
// This function is called in case the exp() operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testExp( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the exp2() operation.
//
// \return void
// \exception std::runtime_error Error in exp2() computation detected.
//
// This function tests the exp2() operation by comparing the results of a vectorized and a
// scalar exp2() operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testExp2( blaze::TrueType )
{
   using blaze::exp2;
   using blaze::loada;
   using blaze::storea;

   test_ = "exp2() operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = exp2( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, exp2( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the exp2() operation.
//
// \return void
//
// This function is called in case the exp2() operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testExp2( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the exp10() operation.
//
// \return void
// \exception std::runtime_error Error in exp10() computation detected.
//
// This function tests the exp10() operation by comparing the results of a vectorized and a
// scalar exp10() operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testExp10( blaze::TrueType )
{
   using blaze::exp10;
   using blaze::loada;
   using blaze::storea;

   test_ = "exp10() operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = pow( T(10), a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, exp10( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the exp10() operation.
//
// \return void
//
// This function is called in case the exp10() operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testExp10( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the natural logarithm operation.
//
// \return void
// \exception std::runtime_error Error in natural logarithm computation detected.
//
// This function tests the natural logarithm operation by comparing the results of a vectorized
// and a scalar exponent operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testLog( blaze::TrueType )
{
   using blaze::log;
   using blaze::loada;
   using blaze::storea;

   test_ = "Natural logarithm operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = log( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, log( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the natural logarithm operation.
//
// \return void
//
// This function is called in case the natural logarithm operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testLog( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the binary logarithm operation.
//
// \return void
// \exception std::runtime_error Error in binary logarithm computation detected.
//
// This function tests the binary logarithm operation by comparing the results of a vectorized
// and a scalar exponent operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testLog2( blaze::TrueType )
{
   using blaze::log2;
   using blaze::loada;
   using blaze::storea;

   test_ = "Binary logarithm operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = log2( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, log2( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the binary logarithm operation.
//
// \return void
//
// This function is called in case the binary logarithm operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testLog2( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the common logarithm operation.
//
// \return void
// \exception std::runtime_error Error in common logarithm computation detected.
//
// This function tests the common logarithm operation by comparing the results of a vectorized
// and a scalar exponent operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testLog10( blaze::TrueType )
{
   using blaze::log10;
   using blaze::loada;
   using blaze::storea;

   test_ = "Common logarithm operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = log10( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, log10( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the common logarithm operation.
//
// \return void
//
// This function is called in case the common logarithm operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testLog10( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the sine operation.
//
// \return void
// \exception std::runtime_error Error in sine computation detected.
//
// This function tests the sine operation by comparing the results of a vectorized and a scalar
// sine operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSin( blaze::TrueType )
{
   using blaze::sin;
   using blaze::loada;
   using blaze::storea;

   test_ = "Sine operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = sin( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, sin( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the sine operation.
//
// \return void
//
// This function is called in case the sine operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSin( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse sine operation.
//
// \return void
// \exception std::runtime_error Error in inverse sine computation detected.
//
// This function tests the inverse sine operation by comparing the results of a vectorized and a
// scalar inverse sine operation. In case any error is detected, a \a std::runtime_error exception
// is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAsin( blaze::TrueType )
{
   using blaze::asin;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse sine operation";

   initialize( T(-1), T(1) );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = asin( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, asin( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse sine operation.
//
// \return void
//
// This function is called in case the inverse sine operation is not available for the given data
// type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAsin( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the hyperbolic sine operation.
//
// \return void
// \exception std::runtime_error Error in hyperbolic sine computation detected.
//
// This function tests the hyperbolic sine operation by comparing the results of a vectorized
// and a scalar hyperbolic sine operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSinh( blaze::TrueType )
{
   using blaze::sinh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Hyperbolic sine operation";

   initialize( T(-1), T(1) );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = sinh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, sinh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the hyperbolic sine operation.
//
// \return void
//
// This function is called in case the hyperbolic sine operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSinh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse hyperbolic sine operation.
//
// \return void
// \exception std::runtime_error Error in inverse hyperbolic sine computation detected.
//
// This function tests the inverse hyperbolic sine operation by comparing the results of
// a vectorized and a scalar inverse hyperbolic sine operation. In case any error is detected, a
// \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAsinh( blaze::TrueType )
{
   using blaze::asinh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse hyperbolic sine operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = asinh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, asinh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse hyperbolic sine operation.
//
// \return void
//
// This function is called in case the inverse hyperbolic sine operation is not available for
// the given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAsinh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the cosine operation.
//
// \return void
// \exception std::runtime_error Error in cosine computation detected.
//
// This function tests the cosine operation by comparing the results of a vectorized and a scalar
// cosine operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCos( blaze::TrueType )
{
   using blaze::cos;
   using blaze::loada;
   using blaze::storea;

   test_ = "Cosine operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = cos( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, cos( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the cosine operation.
//
// \return void
//
// This function is called in case the cosine operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCos( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse cosine operation.
//
// \return void
// \exception std::runtime_error Error in inverse cosine computation detected.
//
// This function tests the inverse cosine operation by comparing the results of a vectorized and
// a scalar inverse cosine operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAcos( blaze::TrueType )
{
   using blaze::acos;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse cosine operation";

   initialize( T(-1), T(1) );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = acos( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, acos( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse cosine operation.
//
// \return void
//
// This function is called in case the inverse cosine operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAcos( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the hyperbolic cosine operation.
//
// \return void
// \exception std::runtime_error Error in hyperbolic cosine computation detected.
//
// This function tests the hyperbolic cosine operation by comparing the results of a vectorized
// and a scalar hyperbolic cosine operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCosh( blaze::TrueType )
{
   using blaze::cosh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Hyperbolic cosine operation";

   initialize( T(-1), T(1) );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = cosh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, cosh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the hyperbolic cosine operation.
//
// \return void
//
// This function is called in case the hyperbolic cosine operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testCosh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse hyperbolic cosine operation.
//
// \return void
// \exception std::runtime_error Error in inverse hyperbolic cosine computation detected.
//
// This function tests the inverse hyperbolic cosine operation by comparing the results of
// a vectorized and a scalar inverse hyperbolic cosine operation. In case any error is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAcosh( blaze::TrueType )
{
   using blaze::acosh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse hyperbolic cosine operation";

   initialize( T(1), T(1000) );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = acosh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, acosh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse hyperbolic cosine operation.
//
// \return void
//
// This function is called in case the inverse hyperbolic cosine operation is not available for
// the given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAcosh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the tangent operation.
//
// \return void
// \exception std::runtime_error Error in tangent computation detected.
//
// This function tests the tangent operation by comparing the results of a vectorized and a scalar
// tangent operation. In case any error is detected, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testTan( blaze::TrueType )
{
   using blaze::tan;
   using blaze::loada;
   using blaze::storea;

   test_ = "Tangent operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = tan( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, tan( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the tangent operation.
//
// \return void
//
// This function is called in case the tangent operation is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testTan( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse tangent operation.
//
// \return void
// \exception std::runtime_error Error in inverse tangent computation detected.
//
// This function tests the inverse tangent operation by comparing the results of a vectorized and
// a scalar inverse tangent operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAtan( blaze::TrueType )
{
   using blaze::atan;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse tangent operation";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = atan( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, atan( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse tangent operation.
//
// \return void
//
// This function is called in case the inverse tangent operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAtan( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the hyperbolic tangent operation.
//
// \return void
// \exception std::runtime_error Error in hyperbolic tangent computation detected.
//
// This function tests the hyperbolic tangent operation by comparing the results of a vectorized
// and a scalar hyperbolic tangent operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testTanh( blaze::TrueType )
{
   using blaze::tanh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Hyperbolic tangent operation";

   initialize( T(-1), T(1) );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = tanh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, tanh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the hyperbolic tangent operation.
//
// \return void
//
// This function is called in case the hyperbolic tangent operation is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testTanh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse hyperbolic tangent operation.
//
// \return void
// \exception std::runtime_error Error in inverse hyperbolic tangent computation detected.
//
// This function tests the inverse hyperbolic tangent operation by comparing the results of
// a vectorized and a scalar inverse hyperbolic tangent operation. In case any error is detected,
// a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAtanh( blaze::TrueType )
{
   using blaze::atanh;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse hyperbolic tangent operation";

   initialize( T(-0.95), T(0.95) );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = atanh( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, atanh( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the inverse tangent operation.
//
// \return void
// \exception std::runtime_error Error in inverse tangent computation detected.
//
// This function tests the inverse tangent operation by comparing the results of a vectorized and
// a scalar inverse tangent operation. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAtan2( blaze::TrueType )
{
   using blaze::atan2;
   using blaze::loada;
   using blaze::storea;

   test_ = "Inverse tangent operation";

   initialize( T(1), T(5) );

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = atan2( a_[i], b_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, atan2( loada( a_+i ), loada( b_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse tangent operation.
//
// \return void
//
// This function is called in case the inverse tangent operation is not available for the given
// data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAtan2( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the inverse hyperbolic tangent operation.
//
// \return void
//
// This function is called in case the inverse hyperbolic tangent operation is not available for
// the given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testAtanh( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the error function (\c erf).
//
// \return void
// \exception std::runtime_error Error in error function computation detected.
//
// This function tests the error function (\c erf) by comparing the results of a vectorized and
// a scalar error function. In case any error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testErf( blaze::TrueType )
{
   using blaze::erf;
   using blaze::loada;
   using blaze::storea;

   test_ = "Error function";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = erf( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, erf( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the error function (\c erf).
//
// \return void
//
// This function is called in case the error function is not available for the given data type
// \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testErf( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the complementary error function (\c erfc).
//
// \return void
// \exception std::runtime_error Error in complementary error function computation detected.
//
// This function tests the complementary error function (\c erfc) by comparing the results of a
// vectorized and a scalar error function. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testErfc( blaze::TrueType )
{
   using blaze::erfc;
   using blaze::loada;
   using blaze::storea;

   test_ = "Complementary error function";

   initialize();

   for( size_t i=0UL; i<N; ++i ) {
      c_[i] = erfc( a_[i] );
   }

   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      storea( d_+i, erfc( loada( a_+i ) ) );
   }

   compare( c_, d_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Skipping the test of the complementary error function (\c erf).
//
// \return void
//
// This function is called in case the complementary error function is not available for the
// given data type \a T.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testErfc( blaze::FalseType )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the addition reduction operation (\c sum).
//
// \return void
// \exception std::runtime_error Error in reduction computation detected.
//
// This function tests the addition reduction operation by comparing the results of a vectorized
// and a scalar reduction. In case any error is detected, a \a std::runtime_error exception is
// thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testSum()
{
   using blaze::loada;
   using blaze::sum;

   test_ = "sum() operation";

   initialize();

   T ssum{};
   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      T tmp{};
      for( size_t j=0UL; j<SIMDSIZE; ++j ) {
         tmp += a_[i+j];
      }
      ssum += tmp;
   }

   T vsum{};
   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      vsum += sum( loada( a_+i ) );
   }

   if( !blaze::equal( ssum, vsum ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Failed reduction operation\n"
          << " Details:\n"
          << "   ssum = " << ssum << "\n"
          << "   vsum = " << vsum << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Testing the multiplication reduction operation (\c prod).
//
// \return void
// \exception std::runtime_error Error in reduction computation detected.
//
// This function tests the multiplication reduction operation by comparing the results of a
// vectorized and a scalar reduction. In case any error is detected, a \a std::runtime_error
// exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::testProd()
{
   using blaze::loada;
   using blaze::prod;

   test_ = "prod() operation";

   initialize();

   T sprod{ 1 };
   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      T tmp{ 1 };
      for( size_t j=0UL; j<SIMDSIZE; ++j ) {
         tmp *= a_[i+j];
      }
      sprod *= tmp;
   }

   T vprod{ 1 };
   for( size_t i=0UL; i<N; i+=SIMDSIZE ) {
      vprod *= prod( loada( a_+i ) );
   }

   if( !blaze::equal( sprod, vprod ) ) {
      std::ostringstream oss;
      oss.precision( 20 );
      oss << " Test : " << test_ << "\n"
          << " Error: Failed reduction operation\n"
          << " Details:\n"
          << "   sprod = " << sprod << "\n"
          << "   vprod = " << vprod << "\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  ERROR DETECTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Comparison of the first 256 elements of the two given arrays.
//
// \param expected The array of expected values.
// \param actual The array of actual values.
// \return void
// \exception std::runtime_error Value mismatch detected.
//
// This function compares the first 256 elements of the two given arrays. In case any value of
// the two arrays differs, a \a std::runtime_error exception is thrown.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::compare( const T* expected, const T* actual ) const
{
   for( size_t i=0UL; i<N; ++i ) {
      if( !blaze::equal( expected[i], actual[i] ) ) {
         std::ostringstream oss;
         oss.precision( 20 );
         oss << " Test : " << test_ << "\n"
             << " Error: Value mismatch detected at index " << i << "\n"
             << " Details:\n"
             << "   expected[" << i << "] = " << expected[i] << "\n"
             << "   actual  [" << i << "] = " << actual[i] << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Initialization of all member arrays.
//
// \return void
//
// This function can be called before each single test case to initialize all arrays with random
// values.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::initialize()
{
   initialize( a_ );
   initialize( b_ );
   initialize( c_ );
   initialize( d_ );
   initialize( e_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of all member arrays.
//
// \param min The smallest possible value for a value.
// \param max The largest possible value for a value.
// \return void
//
// This function can be called before each single test case to initialize all arrays with random
// values in the range \f$ [min,max] \f$.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::initialize( T min, T max )
{
   initialize( a_, min, max );
   initialize( b_, min, max );
   initialize( c_, min, max );
   initialize( d_, min, max );
   initialize( e_, min, max );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of a specific member array.
//
// \param array The array to be initialized.
// \return void
//
// This function can be called before each single test case to initialize the given array with
// random values.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::initialize( T* array )
{
   using blaze::randomize;

   for( size_t i=0UL; i<NN; ++i ) {
      randomize( array[i] );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of a specific member array.
//
// \param array The array to be initialized.
// \param min The smallest possible value for a value.
// \param max The largest possible value for a value.
// \return void
//
// This function can be called before each single test case to initialize the given array with
// random values in the range \f$ [min,max] \f$.
*/
template< typename T >  // Data type of the SIMD test
void OperationTest<T>::initialize( T* array, T min, T max )
{
   using blaze::randomize;

   for( size_t i=0UL; i<NN; ++i ) {
      randomize( array[i], min, max );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Testing the SIMD operations of a specific numeric data type.
//
// \return void
*/
template< typename T >  // Data type of the SIMD test
void runTest()
{
   OperationTest<T>();
}
//*************************************************************************************************




//=================================================================================================
//
//  MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Macro for the execution of an SIMD operation test case.
*/
#define RUN_SIMD_OPERATION_TEST( T ) \
   blazetest::mathtest::simd::runTest<T>()
/*! \endcond */
//*************************************************************************************************

} // namespace simd

} // namespace mathtest

} // namespace blazetest

#endif
