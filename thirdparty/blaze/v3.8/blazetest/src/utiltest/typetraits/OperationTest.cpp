//=================================================================================================
/*!
//  \file src/utiltest/typetraits/OperationTest.cpp
//  \brief Source file for the type traits operation test
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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <blaze/util/constraints/DerivedFrom.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/constraints/SameSize.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/TypeTraits.h>
#include <blazetest/utiltest/typetraits/OperationTest.h>


namespace blazetest {

namespace utiltest {

namespace typetraits {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the OperationTest class test.
//
// \exception std::runtime_error Operation error detected.
*/
OperationTest::OperationTest()
{
   testAddConst();
   testAddCV();
   testAddLValueReference();
   testAddPointer();
   testAddRValueReference();
   testAddVolatile();
   testAll();
   testAny();
   testCommonType();
   testDecay();
   testExtent();
   testGetMember();
   testHasGreaterThan();
   testHasLessThan();
   testHasMember();
   testHasSize();
   testHaveSameSize();
   testIsArithmetic();
   testIsArray();
   testIsBaseOf();
   testIsBoolean();
   testIsBuiltin();
   testIsCharacter();
   testIsClass();
   testIsComplex();
   testIsComplexDouble();
   testIsComplexFloat();
   testIsConst();
   testIsConvertible();
   testIsDouble();
   testIsEmpty();
   testIsEnum();
   testIsFloat();
   testIsFloatingPoint();
   testIsInteger();
   testIsIntegral();
   testIsLong();
   testIsLongDouble();
   testIsLValueReference();
   testIsNumeric();
   testIsObject();
   testIsPod();
   testIsPointer();
   testIsReference();
   testIsRValueReference();
   testIsSame();
   testIsStrictlySame();
   testIsShort();
   testIsSigned();
   testIsUnion();
   testIsUnsigned();
   testIsValid();
   testIsVectorizable();
   testIsVoid();
   testIsVolatile();
   testMakeSigned();
   testMakeUnsigned();
   testRank();
   testRemoveAllExtents();
   testRemoveConst();
   testRemoveCV();
   testRemoveExtent();
   testRemoveLValueReference();
   testRemovePointer();
   testRemoveReference();
   testRemoveRValueReference();
   testRemoveVolatile();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST TYPE TRAITS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the \c AddConst type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c AddConst type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testAddConst()
{
   using blaze::AddConst;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddConst<int>::Type, int const );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddConst<int*>::Type, int* const );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddConst<int&>::Type, int& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddConst<int const>::Type, int const );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddConst<int volatile>::Type, int volatile const );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c AddCV type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c AddCV type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testAddCV()
{
   using blaze::AddCV;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddCV<int>::Type, int const volatile );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddCV<int*>::Type, int* const volatile );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddCV<int&>::Type, int& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddCV<int const>::Type, int const volatile );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddCV<int volatile>::Type, int const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c AddLValueReference type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c AddLValueReference type trait. In case
// an error is detected, a compilation error is created.
*/
void OperationTest::testAddLValueReference()
{
   using blaze::AddLValueReference;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddLValueReference<int>::Type, int& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddLValueReference<int const&>::Type, int const& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddLValueReference<int volatile&&>::Type, int volatile& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddLValueReference<int*>::Type, int*& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddLValueReference<int*&>::Type, int*& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddLValueReference<int*&&>::Type, int*& );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c AddPointer type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c AddPointer type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testAddPointer()
{
   using blaze::AddPointer;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddPointer<int>::Type, int* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddPointer<int const>::Type, int const* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddPointer<int*>::Type, int** );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddPointer<int*&>::Type, int** );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c AddRValueReference type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c AddRValueReference type trait. In case
// an error is detected, a compilation error is created.
*/
void OperationTest::testAddRValueReference()
{
   using blaze::AddRValueReference;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddRValueReference<int>::Type, int&& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddRValueReference<int const&>::Type, int const& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddRValueReference<int volatile&&>::Type, int volatile&& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddRValueReference<int*>::Type, int*&& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddRValueReference<int*&>::Type, int*& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddRValueReference<int*&&>::Type, int*&& );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c AddVolatile type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c AddVolatile type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testAddVolatile()
{
   using blaze::AddVolatile;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddVolatile<int>::Type, int volatile );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddVolatile<int*>::Type, int* volatile );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddVolatile<int&>::Type, int& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddVolatile<int volatile>::Type, int volatile );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( AddVolatile<int const>::Type, int const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c All type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c All type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testAll()
{
   using blaze::All;
   using blaze::IsCharacter;
   using blaze::IsIntegral;
   using blaze::IsPointer;

   const bool value1( All<IsIntegral,int,short,long>::value );
   const bool value2( All<IsIntegral,int,float,double >::value );

   using T1 = All<IsPointer,int*,float*>::Type;
   using T2 = All<IsCharacter,char,signed char,wchar_t>;
   using T3 = All<IsPointer,int*,float&>::Type;
   using T4 = All<IsCharacter,char,signed int,wchar_t>;

   BLAZE_STATIC_ASSERT( value1 == true );
   BLAZE_STATIC_ASSERT( value2 == false );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE   ( T1, blaze::TrueType  );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T2, blaze::TrueType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE   ( T3, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T4, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c Any type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c Any type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testAny()
{
   using blaze::Any;
   using blaze::IsCharacter;
   using blaze::IsIntegral;
   using blaze::IsPointer;

   const bool value1( Any<IsIntegral,int,float>::value );
   const bool value2( Any<IsIntegral,float,double>::value );

   using T1 = Any<IsPointer,int&,float*>::Type;
   using T2 = Any<IsCharacter,float,wchar_t>;
   using T3 = Any<IsPointer,int,float&>::Type;
   using T4 = Any<IsCharacter,int,double>;

   BLAZE_STATIC_ASSERT( value1 == true );
   BLAZE_STATIC_ASSERT( value2 == false );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE   ( T1, blaze::TrueType  );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T2, blaze::TrueType  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE   ( T3, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T4, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c CommonType type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c CommonType type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testCommonType()
{
   using blaze::CommonType;

   using T1 = CommonType<short,int>::Type;
   using T2 = CommonType<const double,int&>::Type;
   using T3 = CommonType<char&, volatile int, const float>::Type;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T1, int    );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T2, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T3, float  );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c Decay type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c Decay type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testDecay()
{
   using blaze::Decay;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( blaze::Decay<int>::Type       , int         );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( blaze::Decay<int&>::Type      , int         );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( blaze::Decay<int&&>::Type     , int         );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( blaze::Decay<const int&>::Type, int         );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( blaze::Decay<int[2]>::Type    , int*        );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( blaze::Decay<int(int)>::Type  , int(*)(int) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c Extent type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c Extent type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testExtent()
{
   using blaze::Extent;

   BLAZE_STATIC_ASSERT( ( Extent< int[4], 0 >::value == 4U ) );
   BLAZE_STATIC_ASSERT( ( Extent< int[2][3][4], 0 >::value == 2U ) );
   BLAZE_STATIC_ASSERT( ( Extent< int[2][3][4], 1 >::value == 3U ) );
   BLAZE_STATIC_ASSERT( ( Extent< int[2][3][4], 2 >::value == 4U ) );
   BLAZE_STATIC_ASSERT( ( Extent< int[][2], 0 >::value == 0U ) );
   BLAZE_STATIC_ASSERT( ( Extent< int[][2], 1 >::value == 2U ) );
   BLAZE_STATIC_ASSERT( ( Extent< int*, 0 >::value == 0U ) );
   BLAZE_STATIC_ASSERT( ( Extent< std::vector<int>, 0 >::value == 0U ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c HasGreaterThan type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c HasGreaterThan type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testHasGreaterThan()
{
   using blaze::HasGreaterThan;

   using T1 = HasGreaterThan< int, int >;
   using T2 = HasGreaterThan< const std::string, std::string >;
   using T3 = HasGreaterThan< volatile int*, int* >;
   using T4 = HasGreaterThan< int, blaze::complex<float> >;
   using T5 = HasGreaterThan< std::string, int >;
   using T6 = HasGreaterThan< int*, std::string* >;

   BLAZE_STATIC_ASSERT( T1::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T2::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T3, blaze::TrueType );
   BLAZE_STATIC_ASSERT( T4::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T5::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T6, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c HasLessThan type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c HasLessThan type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testHasLessThan()
{
   using blaze::HasLessThan;

   using T1 = HasLessThan< int, int >;
   using T2 = HasLessThan< const std::string, std::string >;
   using T3 = HasLessThan< volatile int*, int* >;
   using T4 = HasLessThan< int, blaze::complex<float> >;
   using T5 = HasLessThan< std::string, int >;
   using T6 = HasLessThan< int*, std::string* >;

   BLAZE_STATIC_ASSERT( T1::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T2::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T3, blaze::TrueType );
   BLAZE_STATIC_ASSERT( T4::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T5::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T6, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c GET_MEMBER type trait generation macro.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c GET_MEMBER type trait generation macro.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testGetMember()
{
#if !(defined __INTEL_COMPILER) || ( __INTEL_COMPILER >= 1400 )
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( GetDataType<Type5>::Type, Type5::DataType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( GetDataType<Type6>::Type, Type6::DataType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( GetDataType<Type7>::Type, int             );
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c HAS_MEMBER type trait generation macro.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c HAS_MEMBER type trait generation macro.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testHasMember()
{
   BLAZE_STATIC_ASSERT( HasValue<Type1>::value == true  );
   BLAZE_STATIC_ASSERT( HasValue<Type2>::value == true  );
   BLAZE_STATIC_ASSERT( HasValue<Type3>::value == false );
   BLAZE_STATIC_ASSERT( HasValue<Type4>::value == false );
   BLAZE_STATIC_ASSERT( HasValue<Type5>::value == false );
   BLAZE_STATIC_ASSERT( HasValue<Type6>::value == false );

   BLAZE_STATIC_ASSERT( HasCompute<Type1>::value == false );
   BLAZE_STATIC_ASSERT( HasCompute<Type2>::value == false );
   BLAZE_STATIC_ASSERT( HasCompute<Type3>::value == true  );
   BLAZE_STATIC_ASSERT( HasCompute<Type4>::value == true  );
   BLAZE_STATIC_ASSERT( HasCompute<Type5>::value == false );
   BLAZE_STATIC_ASSERT( HasCompute<Type6>::value == false );

#if !(defined __INTEL_COMPILER) || ( __INTEL_COMPILER >= 1400 )
   BLAZE_STATIC_ASSERT( HasDataType<Type1>::value == false );
   BLAZE_STATIC_ASSERT( HasDataType<Type2>::value == false );
   BLAZE_STATIC_ASSERT( HasDataType<Type3>::value == false );
   BLAZE_STATIC_ASSERT( HasDataType<Type4>::value == false );
   BLAZE_STATIC_ASSERT( HasDataType<Type5>::value == true  );
   BLAZE_STATIC_ASSERT( HasDataType<Type6>::value == true  );
#endif
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c HasSize type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c HasSize type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testHasSize()
{
   using blaze::HasSize;

   using T1 = HasSize<int,4>;
   using T2 = HasSize<float,4>;
   using T3 = HasSize<const double,8>;
   using T4 = HasSize<volatile double,2>;
   using T5 = HasSize<const char,8>;
   using T6 = HasSize<unsigned char,4>;

   BLAZE_STATIC_ASSERT( T1::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T2::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T3, blaze::TrueType );
   BLAZE_STATIC_ASSERT( T4::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T5::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T6, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c HaveSameSize type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c HaveSameSize type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testHaveSameSize()
{
   using blaze::HaveSameSize;

   using T1 = blaze::HaveSameSize<int,unsigned int>;
   using T2 = blaze::HaveSameSize<int,unsigned int>;
   using T3 = blaze::HaveSameSize<int,unsigned int>;
   using T4 = blaze::HaveSameSize<char,wchar_t>;
   using T5 = blaze::HaveSameSize<char,wchar_t>;
   using T6 = blaze::HaveSameSize<char,wchar_t>;

   BLAZE_STATIC_ASSERT( T1::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T2::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T3, blaze::TrueType );
   BLAZE_STATIC_ASSERT( T4::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T5::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T6, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsArithmetic type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsArithmetic type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testIsArithmetic()
{
   using blaze::IsArithmetic;

   BLAZE_STATIC_ASSERT( IsArithmetic<int>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsArithmetic<float const>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsArithmetic<short volatile>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsArithmetic<void>::value == false );
   BLAZE_STATIC_ASSERT( IsArithmetic<int*>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsArithmetic<int&>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsArithmetic<Type7>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsArray type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsArray type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsArray()
{
   using blaze::IsArray;

   BLAZE_STATIC_ASSERT( IsArray< int[3] >::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsArray< const int[] >::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsArray< int[][3] >, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsArray< int >::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsArray< int const* >::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsArray< std::vector<int> >, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsBaseOf type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsBaseOf type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsBaseOf()
{
   using blaze::IsBaseOf;

   class A {};
   class B : public A {};
   class C {};

   using T1 = IsBaseOf<A,B>;
   using T2 = IsBaseOf<A,B>;
   using T3 = IsBaseOf<A,B>;
   using T4 = IsBaseOf<A,C>;
   using T5 = IsBaseOf<B,A>;
   using T6 = IsBaseOf<B,A>;

   BLAZE_STATIC_ASSERT( T1::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T2::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T3, blaze::TrueType );
   BLAZE_STATIC_ASSERT( T4::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T5::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T6, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsBoolean type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsBoolean type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testIsBoolean()
{
   using blaze::IsBoolean;

   BLAZE_STATIC_ASSERT( IsBoolean<bool>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsBoolean<const bool>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsBoolean<const volatile bool>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsBoolean<float>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsBoolean<const int>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsBoolean<volatile short>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsBuiltin type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsBuiltin type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testIsBuiltin()
{
   using blaze::IsBuiltin;

   BLAZE_STATIC_ASSERT( IsBuiltin<void>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsBuiltin<float const>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsBuiltin<short volatile>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsBuiltin<std::string>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsBuiltin<int*>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsBuiltin<int&>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsCharacter type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsCharacter type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testIsCharacter()
{
   using blaze::IsCharacter;

   BLAZE_STATIC_ASSERT( IsCharacter<char>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsCharacter<const unsigned char>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsCharacter<const volatile wchar_t>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsCharacter<unsigned short>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsCharacter<const int>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsCharacter<volatile long>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsClass type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsClass type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testIsClass()
{
   using blaze::IsClass;

   BLAZE_STATIC_ASSERT( IsClass<Type7>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsClass<Type7 const>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsClass<std::string volatile>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsClass<int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsClass<Type7&>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsClass<Type7*>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsComplex type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsComplex type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testIsComplex()
{
   using blaze::complex;
   using blaze::IsComplex;

   BLAZE_STATIC_ASSERT( IsComplex< complex<double> >::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsComplex< const complex<float> >::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsComplex< volatile complex<int> >, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsComplex< float >::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsComplex< const double >::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsComplex< const volatile int >, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsComplexDouble type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsComplexDouble type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testIsComplexDouble()
{
   using blaze::complex;
   using blaze::IsComplexDouble;

   BLAZE_STATIC_ASSERT( IsComplexDouble< complex<double> >::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsComplexDouble< const complex<double> >::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsComplexDouble< volatile complex<double> >, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsComplexDouble< double >::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsComplexDouble< const complex<float> >::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsComplexDouble< const volatile complex<int> >, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsComplexFloat type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsComplexFloat type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testIsComplexFloat()
{
   using blaze::complex;
   using blaze::IsComplexFloat;

   BLAZE_STATIC_ASSERT( IsComplexFloat< complex<float> >::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsComplexFloat< const complex<float> >::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsComplexFloat< volatile complex<float> >, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsComplexFloat< float >::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsComplexFloat< const complex<double> >::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsComplexFloat< const volatile complex<int> >, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsConst type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsConst type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsConst()
{
   using blaze::IsConst;

   BLAZE_STATIC_ASSERT( IsConst<const int>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsConst<const volatile int>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsConst<int* const>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsConst<int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsConst<const int*>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsConst<const int* volatile>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsConvertible type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsConvertible type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testIsConvertible()
{
   using blaze::IsConvertible;

   struct A {};
   struct B : public A {};

   struct C {};
   struct D {
      D( const C& ) {}
   };

   using T1 = IsConvertible<int,unsigned int>;
   using T2 = IsConvertible<float,const double>;
   using T3 = IsConvertible<B,A>;
   using T4 = IsConvertible<B*,A*>;
   using T5 = IsConvertible<C,D>;
   using T6 = IsConvertible<char*,std::string>;
   using T7 = IsConvertible<std::string,char*>;
   using T8 = IsConvertible<A,B>;
   using T9 = IsConvertible<A*,B*>;

   BLAZE_STATIC_ASSERT( T1::value == true );
   BLAZE_STATIC_ASSERT( T2::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T3::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T4::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T5, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T6, blaze::TrueType );
   BLAZE_STATIC_ASSERT( T7::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T8::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T9, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsDouble type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsDouble type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsDouble()
{
   using blaze::IsDouble;

   BLAZE_STATIC_ASSERT( IsDouble<double>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsDouble<const double>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsDouble<const volatile double>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsDouble<float>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsDouble<const int>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsDouble<volatile short>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsEmpty type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsEmpty type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsEmpty()
{
   using blaze::IsEmpty;

   struct A {};
   struct B { int i; };

   BLAZE_STATIC_ASSERT( IsEmpty<A>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsEmpty<A volatile>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsEmpty<A const>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsEmpty<int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsEmpty<std::string>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsEmpty<B>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsEnum type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsEnum type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsEnum()
{
   using blaze::IsEnum;

   enum A {};
   enum B : int {};
   enum class C {};
   class D {};

   BLAZE_STATIC_ASSERT( IsEnum<A>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsEnum<const B>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsEnum<volatile C>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsEnum<int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsEnum<double>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsEnum<D>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsFloat type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsFloat type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsFloat()
{
   using blaze::IsFloat;

   BLAZE_STATIC_ASSERT( IsFloat<float>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsFloat<const float>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsFloat<const volatile float>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsFloat<double>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsFloat<const int>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsFloat<volatile short>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsFloatingPoint type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsFloatingPoint type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testIsFloatingPoint()
{
   using blaze::IsFloatingPoint;

   BLAZE_STATIC_ASSERT( IsFloatingPoint<float>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsFloatingPoint<volatile double>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsFloatingPoint<const long double>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsFloatingPoint<int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsFloatingPoint<const short>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsFloatingPoint<volatile wchar_t>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsInteger type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsInteger type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsInteger()
{
   using blaze::IsInteger;

   BLAZE_STATIC_ASSERT( IsInteger<int>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsInteger<const unsigned int>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsInteger<const volatile signed int>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsInteger<unsigned short>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsInteger<const long>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsInteger<volatile float>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsIntegral type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsIntegral type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsIntegral()
{
   using blaze::IsIntegral;

   BLAZE_STATIC_ASSERT( IsIntegral<int>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsIntegral<const char>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsIntegral<volatile short>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsIntegral<float>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsIntegral<const double>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsIntegral<volatile long double>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsLong type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsLong type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsLong()
{
   using blaze::IsLong;

   BLAZE_STATIC_ASSERT( IsLong<long>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsLong<const unsigned long>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsLong<const volatile signed long>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsLong<unsigned short>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsLong<const int>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsLong<volatile float>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsLongDouble type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsLongDouble type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testIsLongDouble()
{
   using blaze::IsLongDouble;

   BLAZE_STATIC_ASSERT( IsLongDouble<long double>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsLongDouble<const long double>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsLongDouble<const volatile long double>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsLongDouble<float>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsLongDouble<const unsigned int>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsLongDouble<volatile const short>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsLValueReference type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsLValueReference type trait. In case
// an error is detected, a compilation error is created.
*/
void OperationTest::testIsLValueReference()
{
   using blaze::IsLValueReference;

   BLAZE_STATIC_ASSERT( IsLValueReference<int&>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsLValueReference<int (&)(int)>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsLValueReference<const Type1&>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsLValueReference<int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsLValueReference<const Type1&&>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsLValueReference<int (Type7::*)(int)>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsNumeric type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsNumeric type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsNumeric()
{
   using blaze::complex;
   using blaze::IsNumeric;

   BLAZE_STATIC_ASSERT( IsNumeric<int>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsNumeric<const double>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsNumeric<volatile complex<float> >, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsNumeric<void>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsNumeric<bool>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsNumeric<const bool>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsObject type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsObject type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsObject()
{
   using blaze::IsObject;

   BLAZE_STATIC_ASSERT( IsObject<int>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsObject<int*>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsObject<int (*)(void)>, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsObject<int (Type7::*)(void)const>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsObject<int&>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsObject<const void>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsObject<int (double)>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsPod type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsPod type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsPod()
{
   using blaze::IsPod;

   struct A {
      int i_;
      double d_;
   };

   struct B {
      virtual ~B() {}
   };

   struct C {
      std::string s_;
   };

   BLAZE_STATIC_ASSERT( IsPod<int>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsPod<double const>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsPod<A volatile>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsPod< std::vector<int> >::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsPod<B>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsPod<C>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsPointer type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsPointer type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsPointer()
{
   using blaze::IsPointer;

   BLAZE_STATIC_ASSERT( IsPointer<char* const>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsPointer<volatile float*>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsPointer<int (*)(long)>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsPointer<int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsPointer<int Type7::*>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsPointer<int (Type7::*)(long)>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsReference type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsReference type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testIsReference()
{
   using blaze::IsReference;

   BLAZE_STATIC_ASSERT( IsReference<int&>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsReference<int const&>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsReference<int (&)(long)>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsReference<int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsReference<double*>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsReference<int (Type7::*)(long)>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsRValueReference type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsRValueReference type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testIsRValueReference()
{
   using blaze::IsRValueReference;

   BLAZE_STATIC_ASSERT( IsRValueReference<int&&>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsRValueReference<const Type7&&>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsRValueReference<volatile Type7&&>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsRValueReference<int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsRValueReference<const Type7&>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsRValueReference<int (&)(long)>, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsRValueReference<int (Type7::*)(int)>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsSame type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsSame type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsSame()
{
   using blaze::IsSame;

   using T1 = IsSame<int,int>;
   using T2 = IsSame<int,const int>;
   using T3 = IsSame<float,volatile float>;
   using T4 = IsSame<char,wchar_t>;
   using T5 = IsSame<char,volatile float>;
   using T6 = IsSame<int,double>;

   BLAZE_STATIC_ASSERT( T1::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T2::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T3, blaze::TrueType );
   BLAZE_STATIC_ASSERT( T4::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T5::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T6, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsStrictlySame type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsStrictlySame type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsStrictlySame()
{
   using blaze::IsStrictlySame;

   using T1 = IsStrictlySame<int,int>;
   using T2 = IsStrictlySame<const double,const double>;
   using T3 = IsStrictlySame<volatile float,volatile float>;
   using T4 = IsStrictlySame<char,wchar_t>;
   using T5 = IsStrictlySame<int,const int>;
   using T6 = IsStrictlySame<float,volatile float>;

   BLAZE_STATIC_ASSERT( T1::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T2::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T3, blaze::TrueType );
   BLAZE_STATIC_ASSERT( T4::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( T5::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( T6, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsShort type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsShort type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsShort()
{
   using blaze::IsShort;

   BLAZE_STATIC_ASSERT( IsShort<short>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsShort<const unsigned short>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsShort<const volatile signed short>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsShort<unsigned int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsShort<const long>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsShort<volatile float>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsSigned type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsSigned type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsSigned()
{
   using blaze::IsSigned;

   BLAZE_STATIC_ASSERT( IsSigned<short>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsSigned<const int>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsSigned<volatile float>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsSigned<unsigned int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsSigned<const unsigned long>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsSigned<Type7>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsUnion type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsUnion type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsUnion()
{
   using blaze::IsUnion;

   union A;

   BLAZE_STATIC_ASSERT( IsUnion<A>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsUnion<A const>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsUnion<A volatile>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsUnion<int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsUnion<double>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsUnion<std::string>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsUnsigned type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsUnsigned type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsUnsigned()
{
   using blaze::IsUnsigned;

   BLAZE_STATIC_ASSERT( IsUnsigned<unsigned short>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsUnsigned<const unsigned int>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsUnsigned<volatile unsigned long>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsUnsigned<float>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsUnsigned<const volatile int>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsUnsigned<Type7>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsValid type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsValid type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsValid()
{
   using blaze::IsValid;
   using blaze::INVALID_TYPE;

   BLAZE_STATIC_ASSERT( IsValid<int>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsValid<float const>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsValid<double volatile>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsValid<INVALID_TYPE>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsValid<INVALID_TYPE const>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsValid<INVALID_TYPE volatile>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsVectorizable type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsVectorizable type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testIsVectorizable()
{
   using blaze::IsVectorizable;

   BLAZE_STATIC_ASSERT( IsVectorizable< int >::value == bool( BLAZE_SSE2_MODE ) );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsVectorizable< const float >::Type, blaze::BoolConstant< BLAZE_SSE_MODE > );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsVectorizable< volatile double >, blaze::BoolConstant< BLAZE_SSE2_MODE > );
   BLAZE_STATIC_ASSERT( IsVectorizable< void >::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsVectorizable< const bool >::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsVectorizable< volatile Type7 >, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsVoid type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsVoid type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsVoid()
{
   using blaze::IsVoid;

   BLAZE_STATIC_ASSERT( IsVoid<void>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsVoid<const void>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsVoid<const volatile void>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsVoid<int>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsVoid<const char>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsVoid<volatile float>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c IsVolatile type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c IsVolatile type trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsVolatile()
{
   using blaze::IsVolatile;

   BLAZE_STATIC_ASSERT( IsVolatile<volatile int>::value == true );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsVolatile<const volatile int>::Type, blaze::TrueType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsVolatile<int* volatile>, blaze::TrueType );
   BLAZE_STATIC_ASSERT( IsVolatile<volatile int*>::value == false );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( IsVolatile<const int>::Type, blaze::FalseType );
   BLAZE_CONSTRAINT_MUST_BE_DERIVED_FROM( IsVolatile<int>, blaze::FalseType );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c MakeSigned type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c MakeSigned type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testMakeSigned()
{
   using blaze::MakeSigned;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<signed char>::Type   , signed char );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<unsigned char>::Type , signed char );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<short>::Type         , short       );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<unsigned short>::Type, short       );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<int>::Type           , int         );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<unsigned int>::Type  , int         );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<long>::Type          , long        );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<unsigned long>::Type , long        );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<const int>::Type         , const int          );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<volatile int>::Type      , volatile int       );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeSigned<const volatile int>::Type, const volatile int );

   BLAZE_CONSTRAINT_MUST_HAVE_SAME_SIZE( MakeSigned<wchar_t>::Type, wchar_t );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c MakeUnsigned type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c MakeUnsigned type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testMakeUnsigned()
{
   using blaze::MakeUnsigned;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<signed char>::Type   , unsigned char  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<unsigned char>::Type , unsigned char  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<short>::Type         , unsigned short );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<unsigned short>::Type, unsigned short );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<int>::Type           , unsigned int   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<unsigned int>::Type  , unsigned int   );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<long>::Type          , unsigned long  );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<unsigned long>::Type , unsigned long  );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<const int>::Type         , const unsigned int          );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<volatile int>::Type      , volatile unsigned int       );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeUnsigned<const volatile int>::Type, const volatile unsigned int );

   BLAZE_CONSTRAINT_MUST_HAVE_SAME_SIZE( MakeUnsigned<wchar_t>::Type, wchar_t );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c Rank type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c Rank type trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testRank()
{
   using blaze::Rank;

   BLAZE_STATIC_ASSERT( Rank< int[] >::value == 1UL );
   BLAZE_STATIC_ASSERT( Rank< int[3] >::value == 1UL );
   BLAZE_STATIC_ASSERT( Rank< const int[2][3][4] >::value == 3UL );
   BLAZE_STATIC_ASSERT( Rank< int[][3] >::value == 2UL );
   BLAZE_STATIC_ASSERT( Rank< int const* >::value == 0UL );
   BLAZE_STATIC_ASSERT( Rank< std::vector<int> >::value == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c RemoveAllExtents type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c RemoveAllExtents type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testRemoveAllExtents()
{
   using blaze::RemoveAllExtents;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveAllExtents<int>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveAllExtents<int const[2]>::Type, const int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveAllExtents<int[2][4]>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveAllExtents<int[][2]>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveAllExtents<int[2][3][4]>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveAllExtents<int const*>::Type, int const* );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c RemoveConst type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c RemoveConst type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testRemoveConst()
{
   using blaze::RemoveConst;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveConst<short>::Type, short );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveConst<const double>::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveConst<const volatile int>::Type, volatile int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveConst<int const*>::Type, const int* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveConst<int const* const>::Type, const int* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveConst<int const&>::Type, const int& );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c RemoveCV type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c RemoveCV type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testRemoveCV()
{
   using blaze::RemoveCV;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveCV<short>::Type, short );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveCV<const double>::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveCV<volatile float>::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveCV<const volatile int>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveCV<int const*>::Type, int const* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveCV<int const* const>::Type, int const* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveCV<int const&>::Type, int const& );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c RemoveExtent type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c RemoveExtent type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testRemoveExtent()
{
   using blaze::RemoveExtent;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveExtent<int>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveExtent<int const[2]>::Type, const int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveExtent<int[2][4]>::Type, int[4] );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveExtent<int[][2]>::Type, int[2] );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveExtent<int const*>::Type, const int* );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c RemoveLValueReference type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c RemoveLValueReference type trait. In
// case an error is detected, a compilation error is created.
*/
void OperationTest::testRemoveLValueReference()
{
   using blaze::RemoveLValueReference;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveLValueReference<int>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveLValueReference<const int&>::Type, const int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveLValueReference<volatile int&&>::Type, volatile int&& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveLValueReference<int*>::Type, int* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveLValueReference<int*&>::Type, int* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveLValueReference<int*&&>::Type, int*&& );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c RemovePointer type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c RemovePointer type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testRemovePointer()
{
   using blaze::RemovePointer;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemovePointer<int>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemovePointer<const int*>::Type, const int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemovePointer<volatile int**>::Type, volatile int* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemovePointer<int&>::Type, int& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemovePointer<int*&>::Type, int*& );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c RemoveReference type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c RemoveReference type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testRemoveReference()
{
   using blaze::RemoveReference;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveReference<int>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveReference<const int&>::Type, const int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveReference<volatile int&&>::Type, volatile int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveReference<int*>::Type, int* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveReference<int*&>::Type, int* );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c RemoveRValueReference type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c RemoveRValueReference type trait. In
// case an error is detected, a compilation error is created.
*/
void OperationTest::testRemoveRValueReference()
{
   using blaze::RemoveRValueReference;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveRValueReference<int>::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveRValueReference<const int&>::Type, const int& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveRValueReference<volatile int&&>::Type, volatile int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveRValueReference<int*>::Type, int* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveRValueReference<int*&>::Type, int*& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveRValueReference<int*&&>::Type, int* );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c RemoveVolatile type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the \c RemoveVolatile type trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testRemoveVolatile()
{
   using blaze::RemoveVolatile;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveVolatile<short>::Type, short );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveVolatile<volatile double>::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveVolatile<const volatile int>::Type, const int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveVolatile<int volatile*>::Type, int volatile* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveVolatile<int volatile* volatile>::Type, int volatile* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( RemoveVolatile<int volatile&>::Type, int volatile& );
}
//*************************************************************************************************

} // namespace typetraits

} // namespace utiltest

} // namespace blazetest




//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running type traits operation test..." << std::endl;

   try
   {
      RUN_TYPETRAITS_OPERATION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during type traits operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
