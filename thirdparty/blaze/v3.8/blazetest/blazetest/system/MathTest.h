//=================================================================================================
/*!
//  \file blazetest/system/MathTest.h
//  \brief General settings for the math tests of the blaze test suite
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

#ifndef _BLAZETEST_SYSTEM_MATHTEST_H_
#define _BLAZETEST_SYSTEM_MATHTEST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/typetraits/HasAdd.h>
#include <blaze/math/typetraits/HasDiv.h>
#include <blaze/math/typetraits/HasMax.h>
#include <blaze/math/typetraits/HasMin.h>
#include <blaze/math/typetraits/HasMult.h>
#include <blaze/math/typetraits/HasSub.h>
#include <blaze/math/typetraits/IsBLASCompatible.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/util/constraints/Valid.h>
#include <blaze/util/StaticAssert.h>
#include <blazetest/system/Types.h>


namespace blazetest {

namespace mathtest {

//=================================================================================================
//
//  USING DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
using blaze::rowVector;
using blaze::columnVector;
using blaze::rowMajor;
using blaze::columnMajor;

using blaze::CompressedMatrix;
using blaze::CompressedVector;
using blaze::DynamicMatrix;
using blaze::DynamicVector;
using blaze::StaticMatrix;
using blaze::StaticVector;
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GENERAL CONFIGURATION
//
//=================================================================================================

#include <blazetest/config/MathTest.h>




//=================================================================================================
//
//  BASIC DATA TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief First data type for the math tests.
//
// This type represents the first basic element type used for vector and matrix operation tests.
*/
using TypeA = BLAZETEST_MATHTEST_TYPE_A;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Second data type for the math tests.
//
// This type represents the second basic element type used for vector and matrix operation tests.
*/
using TypeB = BLAZETEST_MATHTEST_TYPE_B;
//*************************************************************************************************




//=================================================================================================
//
//  DERIVED DATA TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The scalar element type of TypeA.
//
// This type represents the underlying scalar element type of the specified TypeA. It is used
// for vector and matrix type that only support scalar data types.
*/
using ScalarA = blaze::UnderlyingScalar_t<TypeA>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The scalar element type of TypeB.
//
// This type represents the underlying scalar element type of the specified TypeB. It is used
// for vector and matrix type that only support scalar data types.
*/
using ScalarB = blaze::UnderlyingScalar_t<TypeB>;
//*************************************************************************************************




//=================================================================================================
//
//  COMPILE TIME CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace {

BLAZE_STATIC_ASSERT( blaze::IsScalar<TypeA>::value || blaze::IsVector<TypeA>::value || blaze::IsMatrix<TypeA>::value );
BLAZE_STATIC_ASSERT( blaze::IsScalar<TypeB>::value || blaze::IsVector<TypeB>::value || blaze::IsMatrix<TypeB>::value );
BLAZE_STATIC_ASSERT( blaze::IsScalar<ScalarA>::value );
BLAZE_STATIC_ASSERT( blaze::IsScalar<ScalarB>::value );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_ADDITION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_ADDITION > 2 ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_ADDITION || ( blaze::HasAdd<TypeA,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_ADDITION || ( blaze::HasAdd<TypeA,TypeB>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_ADDITION || ( blaze::HasAdd<TypeB,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_ADDITION || ( blaze::HasAdd<TypeB,TypeB>::value ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SUBTRACTION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SUBTRACTION > 2 ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_SUBTRACTION || ( blaze::HasSub<TypeA,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_SUBTRACTION || ( blaze::HasSub<TypeA,TypeB>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_SUBTRACTION || ( blaze::HasSub<TypeB,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_SUBTRACTION || ( blaze::HasSub<TypeB,TypeB>::value ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_MULTIPLICATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_MULTIPLICATION > 2 ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MULTIPLICATION || ( blaze::HasMult<TypeA,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MULTIPLICATION || ( blaze::HasMult<TypeA,TypeB>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MULTIPLICATION || ( blaze::HasMult<TypeB,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MULTIPLICATION || ( blaze::HasMult<TypeB,TypeB>::value ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DIVISION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DIVISION > 2 ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_DIVISION || ( blaze::HasDiv<TypeA,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_DIVISION || ( blaze::HasDiv<TypeA,TypeB>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_DIVISION || ( blaze::HasDiv<TypeB,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_DIVISION || ( blaze::HasDiv<TypeB,TypeB>::value ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_MINIMUM < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_MINIMUM > 2 ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MINIMUM || ( blaze::HasMin<TypeA,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MINIMUM || ( blaze::HasMin<TypeA,TypeB>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MINIMUM || ( blaze::HasMin<TypeB,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MINIMUM || ( blaze::HasMin<TypeB,TypeB>::value ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_MAXIMUM < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_MAXIMUM > 2 ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MAXIMUM || ( blaze::HasMax<TypeA,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MAXIMUM || ( blaze::HasMax<TypeA,TypeB>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MAXIMUM || ( blaze::HasMax<TypeB,TypeA>::value ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_MAXIMUM || ( blaze::HasMax<TypeB,TypeB>::value ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_BASIC_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SCALED_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SCALED_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_TRANS_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_TRANS_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_ABS_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_ABS_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_CONJ_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_CONJ_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_REAL_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_REAL_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_IMAG_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_IMAG_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_INV_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_INV_OPERATION > 2 ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_INV_OPERATION || blaze::IsBLASCompatible<TypeA>::value );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_INV_OPERATION || blaze::IsBLASCompatible<TypeB>::value );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_EVAL_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_EVAL_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SERIAL_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SERIAL_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_NOALIAS_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_NOALIAS_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_NOSIMD_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_NOSIMD_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION > 2 ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION || blaze::IsScalar<TypeA>::value );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION || blaze::IsScalar<TypeB>::value );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION > 2 ) );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION || blaze::IsScalar<TypeA>::value );
BLAZE_STATIC_ASSERT( !BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION || blaze::IsScalar<TypeB>::value );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_ROW_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_ROW_OPERATION > 2 ) );

BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION < 0 ) );
BLAZE_STATIC_ASSERT( !( BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION > 2 ) );

}
/*! \endcond */
//*************************************************************************************************

} // namespace mathtest

} // namespace blazetest

#endif
