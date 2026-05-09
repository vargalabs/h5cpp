//=================================================================================================
/*!
//  \file src/mathtest/typetraits/OperationTest.cpp
//  \brief Source file for the mathematical type traits operation test
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
#include <vector>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/Commutative.h>
#include <blaze/math/constraints/CompositeType.h>
#include <blaze/math/constraints/CUDAAssignable.h>
#include <blaze/math/constraints/Diagonal.h>
#include <blaze/math/constraints/Hermitian.h>
#include <blaze/math/constraints/Identity.h>
#include <blaze/math/constraints/Invertible.h>
#include <blaze/math/constraints/Lower.h>
#include <blaze/math/constraints/Matrix.h>
#include <blaze/math/constraints/PaddingEnabled.h>
#include <blaze/math/constraints/ResultType.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/SIMDEnabled.h>
#include <blaze/math/constraints/Static.h>
#include <blaze/math/constraints/StrictlyLower.h>
#include <blaze/math/constraints/StrictlyUpper.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/Uniform.h>
#include <blaze/math/constraints/UniLower.h>
#include <blaze/math/constraints/UniUpper.h>
#include <blaze/math/constraints/Upper.h>
#include <blaze/math/constraints/Vector.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/IdentityMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/StrictlyLowerMatrix.h>
#include <blaze/math/StrictlyUpperMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/typetraits/DynamicAllocator.h>
#include <blaze/math/typetraits/GetAllocator.h>
#include <blaze/math/typetraits/HasCompositeType.h>
#include <blaze/math/typetraits/HasResultType.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsCommutative.h>
#include <blaze/math/typetraits/IsCUDAAssignable.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsInvertible.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsPaddingEnabled.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsSIMDEnabled.h>
#include <blaze/math/typetraits/IsStatic.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/MakeComplex.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/math/typetraits/UnderlyingNumeric.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/math/ZeroMatrix.h>
#include <blaze/util/constraints/SameType.h>
#include <blazetest/mathtest/typetraits/OperationTest.h>


namespace blazetest {

namespace mathtest {

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
   testGetAllocator();
   testHasCompositeType();
   testHasResultType();
   testIsColumnMajorMatrix();
   testIsColumnVector();
   testIsCommutative();
   testIsCUDAAssignable();
   testIsDiagonal();
   testIsHermitian();
   testIsIdentity();
   testIsInvertible();
   testIsLower();
   testIsMatrix();
   testIsPaddingEnabled();
   testIsRowVector();
   testIsScalar();
   testIsSIMDEnabled();
   testIsStatic();
   testIsStrictlyLower();
   testIsStrictlyUpper();
   testIsSymmetric();
   testIsUniform();
   testIsUniLower();
   testIsUniUpper();
   testIsUpper();
   testIsVector();
   testIsZero();
   testMakeComplex();
   testRemoveAdaptor();
   testUnderlyingBuiltin();
   testUnderlyingElement();
   testUnderlyingScalar();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST TYPE TRAITS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the mathematical 'DynamicAllocator' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'DynamicAllocator' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testDynamicAllocator()
{
   using blaze::DynamicAllocator;
   using blaze::AlignedAllocator;
   using blaze::NullAllocator;

   using Alloc1 = AlignedAllocator<int>;
   using Alloc2 = NullAllocator<double>;

   using Result1 = DynamicAllocator<Alloc1>::Type<A>;
   using Result2 = DynamicAllocator<Alloc2>::Type<A>;
   using Result3 = DynamicAllocator<Alloc1,Alloc1>::Type<A>;
   using Result4 = DynamicAllocator<Alloc1,Alloc2>::Type<A>;
   using Result5 = DynamicAllocator<Alloc2,Alloc1>::Type<A>;
   using Result6 = DynamicAllocator<Alloc2,Alloc2>::Type<A>;

   using Expected = AlignedAllocator<A>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Result1, Expected );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Result2, Expected );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Result3, Expected );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Result4, Expected );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Result5, Expected );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Result6, Expected );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'GetAllocator' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'GetAllocator' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testGetAllocator()
{
   using blaze::DynamicVector;
   using blaze::DynamicMatrix;
   using blaze::StaticVector;
   using blaze::CompressedMatrix;
   using blaze::GetAllocator;
   using blaze::AlignedAllocator;
   using blaze::NullAllocator;

   using Source1 = DynamicVector<int>;
   using Source2 = const DynamicVector<double>;
   using Source3 = volatile DynamicMatrix<int>;
   using Source4 = int;
   using Source5 = const StaticVector<float,3UL>;
   using Source6 = volatile CompressedMatrix<int>;

   using Result1 = AlignedAllocator<int>;
   using Result2 = AlignedAllocator<double>;
   using Result3 = AlignedAllocator<int>;
   using Result4 = NullAllocator<int>;
   using Result5 = NullAllocator<float>;
   using Result6 = NullAllocator<int>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( GetAllocator<Source1>::Type, Result1 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( GetAllocator<Source2>::Type, Result2 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( GetAllocator<Source3>::Type, Result3 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( GetAllocator<Source4>::Type, Result4 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( GetAllocator<Source5>::Type, Result5 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( GetAllocator<Source6>::Type, Result6 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'HasCompositeType' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'HasCompositeType' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testHasCompositeType()
{
   using Type1 = int;
   using Type2 = blaze::complex<float>;
   using Type3 = blaze::DynamicVector<int>;
   using Type4 = blaze::CompressedVector<int>;
   using Type5 = blaze::DynamicMatrix<int>;
   using Type6 = blaze::CompressedMatrix<int>;

   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type3                 );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type3 const           );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type4                 );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type4 const           );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type5                 );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type5 const           );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type6                 );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type6 const           );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_HAVE_COMPOSITE_TYPE    ( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_COMPOSITE_TYPE( Type6* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'HasResultType' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'HasResultType' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testHasResultType()
{
   using Type1 = int;
   using Type2 = blaze::complex<float>;
   using Type3 = blaze::DynamicVector<int>;
   using Type4 = blaze::CompressedVector<int>;
   using Type5 = blaze::DynamicMatrix<int>;
   using Type6 = blaze::CompressedMatrix<int>;

   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type3                 );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type3 const           );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type4                 );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type4 const           );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type5                 );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type5 const           );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type6                 );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type6 const           );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_HAVE_RESULT_TYPE    ( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_HAVE_RESULT_TYPE( Type6* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsColumnMajorMatrix' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsColumnMajorMatrix' type
// trait. In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsColumnMajorMatrix()
{
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::columnVector;
   using blaze::rowVector;

   using Type1 = blaze::DynamicMatrix<double,columnMajor>;
   using Type2 = blaze::CompressedMatrix<int,columnMajor>;
   using Type3 = int;
   using Type4 = blaze::complex<float>;
   using Type5 = blaze::DynamicVector<double,columnVector>;
   using Type6 = blaze::CompressedVector<int,rowVector>;
   using Type7 = blaze::DynamicMatrix<double,rowMajor>;
   using Type8 = blaze::CompressedMatrix<int,rowMajor>;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE    ( Type1                 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE    ( Type1 const           );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE    ( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE    ( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE    ( Type2                 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE    ( Type2 const           );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE    ( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE    ( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type5                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type5 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type6                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type6 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type6* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type7                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type7 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type7 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type7 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type7&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type7*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type7* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type7* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type7* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type8                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type8 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type8 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type8 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type8&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type8*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type8* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type8* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_MAJOR_MATRIX_TYPE( Type8* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsColumnVector' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsColumnVector' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsColumnVector()
{
   using blaze::columnVector;
   using blaze::rowVector;
   using blaze::rowMajor;
   using blaze::columnMajor;

   using Type1 = blaze::DynamicVector<double,columnVector>;
   using Type2 = blaze::CompressedVector<int,columnVector>;
   using Type3 = int;
   using Type4 = blaze::complex<float>;
   using Type5 = blaze::DynamicVector<double,rowVector>;
   using Type6 = blaze::CompressedVector<int,rowVector>;
   using Type7 = blaze::DynamicMatrix<double,rowMajor>;
   using Type8 = blaze::CompressedMatrix<int,columnMajor>;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE    ( Type1                 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE    ( Type1 const           );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE    ( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE    ( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE    ( Type2                 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE    ( Type2 const           );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE    ( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE    ( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type5                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type5 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type6                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type6 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type6* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type7                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type7 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type7 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type7 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type7&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type7*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type7* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type7* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type7* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type8                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type8 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type8 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type8 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type8&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type8*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type8* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type8* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMN_VECTOR_TYPE( Type8* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsCommutative' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsCommutative' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsCommutative()
{
   using VT = blaze::StaticVector<int,3UL>;
   using MT = blaze::StaticMatrix<int,3UL,3UL>;

   using Type1 = double;
   using Type2 = blaze::complex<double>;
   using Type3 = blaze::DynamicVector<int>;
   using Type4 = blaze::DynamicVector<VT>;
   using Type5 = blaze::DynamicVector<MT>;
   using Type6 = blaze::DynamicMatrix<int>;
   using Type7 = blaze::DynamicMatrix<VT>;
   using Type8 = blaze::DynamicMatrix<MT>;

   BLAZE_CONSTRAINT_MUST_BE_COMMUTATIVE_TYPES    ( Type1, Type2 );
   BLAZE_CONSTRAINT_MUST_BE_COMMUTATIVE_TYPES    ( Type3, Type3 );
   BLAZE_CONSTRAINT_MUST_BE_COMMUTATIVE_TYPES    ( Type7, Type7 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMMUTATIVE_TYPES( Type6, Type3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMMUTATIVE_TYPES( Type5, Type4 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMMUTATIVE_TYPES( Type7, Type8 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsCUDAAssignable' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsCUDAAssignable' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsCUDAAssignable()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_CUDA_ASSIGNABLE( A );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CUDA_ASSIGNABLE( I );
   BLAZE_CONSTRAINT_MUST_BE_CUDA_ASSIGNABLE    ( J );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsDiagonal' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsDiagonal' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsDiagonal()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::LowerMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;
   using blaze::IdentityMatrix;

   using Type1 = SymmetricMatrix< DynamicMatrix<int> >;
   using Type2 = LowerMatrix< DynamicMatrix<int> >;
   using Type3 = UpperMatrix< DynamicMatrix<int> >;
   using Type4 = DiagonalMatrix< DynamicMatrix<int> >;
   using Type5 = IdentityMatrix<int>;

   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE    ( Type4                 );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE    ( Type4 const           );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE    ( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE    ( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE    ( Type5                 );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE    ( Type5 const           );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE    ( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_DIAGONAL_MATRIX_TYPE    ( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( Type5* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsHermitian' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsHermitian' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsHermitian()
{
   using blaze::complex;
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::HermitianMatrix;
   using blaze::DiagonalMatrix;
   using blaze::IdentityMatrix;

   using Type1 = DynamicMatrix<int>;
   using Type2 = SymmetricMatrix< DynamicMatrix<int> >;
   using Type3 = SymmetricMatrix< DynamicMatrix< complex<int> > >;
   using Type4 = HermitianMatrix< DynamicMatrix<int> >;
   using Type5 = HermitianMatrix< DynamicMatrix< complex<int> > >;
   using Type6 = DiagonalMatrix< DynamicMatrix<int> >;
   using Type7 = DiagonalMatrix< DynamicMatrix< complex<int> > >;
   using Type8 = IdentityMatrix<int>;
   using Type9 = IdentityMatrix< complex<int> >;

   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type2                 );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type2 const           );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type4                 );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type4 const           );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type5                 );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type5 const           );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type6                 );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type6 const           );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type6* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type7                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type7 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type7 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type7 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type7&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type7*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type7* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type7* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type7* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type8                 );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type8 const           );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type8 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type8 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type8&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type8*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type8* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type8* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type8* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type9                 );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type9 const           );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type9 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_HERMITIAN_MATRIX_TYPE    ( Type9 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type9&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type9*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type9* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type9* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( Type9* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsIdentity' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsIdentity' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsIdentity()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::LowerMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;
   using blaze::IdentityMatrix;

   using Type1 = SymmetricMatrix< DynamicMatrix<int> >;
   using Type2 = LowerMatrix< DynamicMatrix<int> >;
   using Type3 = UpperMatrix< DynamicMatrix<int> >;
   using Type4 = DiagonalMatrix< DynamicMatrix<int> >;
   using Type5 = IdentityMatrix<int>;

   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_IDENTITY_MATRIX_TYPE    ( Type5                 );
   BLAZE_CONSTRAINT_MUST_BE_IDENTITY_MATRIX_TYPE    ( Type5 const           );
   BLAZE_CONSTRAINT_MUST_BE_IDENTITY_MATRIX_TYPE    ( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_IDENTITY_MATRIX_TYPE    ( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_IDENTITY_MATRIX_TYPE( Type5* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsInvertible' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsInvertible' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsInvertible()
{
   using Type1  = float;
   using Type2  = double;
   using Type3  = long double;
   using Type4  = blaze::complex<float>;
   using Type5  = blaze::complex<double>;
   using Type6  = blaze::complex<long double>;
   using Type7  = blaze::DynamicMatrix<double>;
   using Type8  = int;
   using Type9  = blaze::complex<int>;
   using Type10 = blaze::DynamicMatrix<int>;

   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type1                  );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type1 const            );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type1 volatile         );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type1 const volatile   );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type1&                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type1*                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type1* const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type1* volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type1* const volatile  );

   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type2                  );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type2 const            );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type2 volatile         );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type2 const volatile   );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type2&                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type2*                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type2* const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type2* volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type2* const volatile  );

   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type3                  );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type3 const            );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type3 volatile         );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type3 const volatile   );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type3&                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type3*                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type3* const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type3* volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type3* const volatile  );

   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type4                  );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type4 const            );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type4 volatile         );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type4 const volatile   );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type4&                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type4*                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type4* const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type4* volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type4* const volatile  );

   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type5                  );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type5 const            );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type5 volatile         );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type5 const volatile   );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type5&                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type5*                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type5* const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type5* volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type5* const volatile  );

   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type6                  );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type6 const            );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type6 volatile         );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type6 const volatile   );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type6&                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type6*                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type6* const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type6* volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type6* const volatile  );

   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type7                  );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type7 const            );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type7 volatile         );
   BLAZE_CONSTRAINT_MUST_BE_INVERTIBLE_TYPE    ( Type7 const volatile   );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type7&                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type7*                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type7* const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type7* volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type7* const volatile  );

   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type8                  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type8 const            );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type8 volatile         );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type8 const volatile   );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type8&                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type8*                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type8* const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type8* volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type8* const volatile  );

   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type9                  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type9 const            );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type9 volatile         );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type9 const volatile   );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type9&                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type9*                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type9* const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type9* volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type9* const volatile  );

   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type10                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type10 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type10 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type10 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type10&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type10*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type10* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type10* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_INVERTIBLE_TYPE( Type10* const volatile );

}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsLower' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsLower' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsLower()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::LowerMatrix;
   using blaze::UniLowerMatrix;
   using blaze::StrictlyLowerMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;
   using blaze::IdentityMatrix;

   using Type1  = DynamicMatrix<int>;
   using Type2  = SymmetricMatrix< DynamicMatrix<int> >;
   using Type3  = LowerMatrix< DynamicMatrix<int> >;
   using Type4  = UniLowerMatrix< DynamicMatrix<int> >;
   using Type5  = StrictlyLowerMatrix< DynamicMatrix<int> >;
   using Type6  = UpperMatrix< DynamicMatrix<int> >;
   using Type7  = DiagonalMatrix< DynamicMatrix<int> >;
   using Type8  = IdentityMatrix<int>;

   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type3                 );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type3 const           );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type4                 );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type4 const           );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type5                 );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type5 const           );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type6                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type6 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type6* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type7                 );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type7 const           );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type7 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type7 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type7&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type7*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type7* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type7* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type7* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type8                 );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type8 const           );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type8 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_LOWER_MATRIX_TYPE    ( Type8 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type8&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type8*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type8* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type8* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( Type8* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsMatrix' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsMatrix' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsMatrix()
{
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::columnVector;
   using blaze::rowVector;

   using Type1 = blaze::DynamicMatrix<double,rowMajor>;
   using Type2 = blaze::CompressedMatrix<int,columnMajor>;
   using Type3 = int;
   using Type4 = blaze::complex<float>;
   using Type5 = blaze::DynamicVector<double,columnVector>;
   using Type6 = blaze::CompressedVector<int,rowVector>;

   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE    ( Type1                 );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE    ( Type1 const           );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE    ( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE    ( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE    ( Type2                 );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE    ( Type2 const           );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE    ( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE    ( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type5                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type5 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type6                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type6 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_MATRIX_TYPE( Type6* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsPaddingEnabled' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsPaddingEnabled' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsPaddingEnabled()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_PADDING_ENABLED( A );
   BLAZE_CONSTRAINT_MUST_NOT_BE_PADDING_ENABLED( E );
   BLAZE_CONSTRAINT_MUST_BE_PADDING_ENABLED    ( F );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsRowMajorMatrix' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsRowMajorMatrix' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsRowMajorMatrix()
{
   using blaze::rowMajor;
   using blaze::columnMajor;
   using blaze::columnVector;
   using blaze::rowVector;

   using Type1 = blaze::DynamicMatrix<double,rowMajor>;
   using Type2 = blaze::CompressedMatrix<int,rowMajor>;
   using Type3 = int;
   using Type4 = blaze::complex<float>;
   using Type5 = blaze::DynamicVector<double,columnVector>;
   using Type6 = blaze::CompressedVector<int,rowVector>;
   using Type7 = blaze::DynamicMatrix<double,columnMajor>;
   using Type8 = blaze::CompressedMatrix<int,columnMajor>;

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( Type1                 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( Type1 const           );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( Type2                 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( Type2 const           );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type5                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type5 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type6                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type6 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type6* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type7                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type7 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type7 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type7 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type7&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type7*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type7* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type7* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type7* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type8                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type8 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type8 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type8 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type8&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type8*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type8* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type8* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_MAJOR_MATRIX_TYPE( Type8* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsRowVector' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsRowVector' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsRowVector()
{
   using blaze::columnVector;
   using blaze::rowVector;
   using blaze::rowMajor;
   using blaze::columnMajor;

   using Type1 = blaze::DynamicVector<double,rowVector>;
   using Type2 = blaze::CompressedVector<int,rowVector>;
   using Type3 = int;
   using Type4 = blaze::complex<float>;
   using Type5 = blaze::DynamicVector<double,columnVector>;
   using Type6 = blaze::CompressedVector<int,columnVector>;
   using Type7 = blaze::DynamicMatrix<double,rowMajor>;
   using Type8 = blaze::CompressedMatrix<int,columnMajor>;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( Type1                 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( Type1 const           );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( Type2                 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( Type2 const           );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type5                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type5 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type6                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type6 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type6* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type7                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type7 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type7 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type7 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type7&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type7*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type7* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type7* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type7* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type8                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type8 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type8 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type8 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type8&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type8*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type8* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type8* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROW_VECTOR_TYPE( Type8* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsScalar' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsScalar' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsScalar()
{
   using blaze::columnVector;
   using blaze::rowVector;
   using blaze::rowMajor;
   using blaze::columnMajor;

   using Type1 = int;
   using Type2 = blaze::complex<float>;
   using Type3 = blaze::DynamicVector<float,columnVector>;
   using Type4 = blaze::CompressedVector<int,rowVector>;
   using Type5 = blaze::DynamicMatrix<double,rowMajor>;
   using Type6 = blaze::CompressedMatrix<int,columnMajor>;

   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE    ( Type1                 );
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE    ( Type1 const           );
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE    ( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE    ( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE    ( Type2                 );
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE    ( Type2 const           );
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE    ( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE    ( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type5                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type5 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type6                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type6 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SCALAR_TYPE( Type6* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsSIMDEnabled' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsSIMDEnabled' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsSIMDEnabled()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SIMD_ENABLED( A );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SIMD_ENABLED( G );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_ENABLED    ( H );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsStatic' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsStatic' type trait. In
// case an error is detected, a compilation error is created.
*/
void OperationTest::testIsStatic()
{
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::StaticMatrix;
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;

   using Type1 = StaticVector<int,3UL>;
   using Type2 = DynamicVector<int>;
   using Type3 = StaticMatrix<int,3UL,4UL>;
   using Type4 = DynamicMatrix<int>;
   using Type5 = SymmetricMatrix< StaticMatrix<int,3UL,3UL> >;
   using Type6 = SymmetricMatrix< DynamicMatrix<int> >;

   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type1                 );
   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type1 const           );
   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type3                 );
   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type3 const           );
   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type5                 );
   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type5 const           );
   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_STATIC_TYPE    ( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type6                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type6 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STATIC_TYPE( Type6* const volatile );
}
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Test of the mathematical 'IsStrictlyLower' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsStrictlyLower' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsStrictlyLower()
{
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;
   using blaze::StrictlyLowerMatrix;
   using blaze::StrictlyUpperMatrix;
   using blaze::DiagonalMatrix;

   using Type1 = DynamicMatrix<int>;
   using Type2 = LowerMatrix< DynamicMatrix<int> >;
   using Type3 = StrictlyLowerMatrix< DynamicMatrix<int> >;
   using Type4 = StrictlyUpperMatrix< DynamicMatrix<int> >;
   using Type5 = DiagonalMatrix< DynamicMatrix<int> >;

   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_LOWER_MATRIX_TYPE    ( Type3                 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_LOWER_MATRIX_TYPE    ( Type3 const           );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_LOWER_MATRIX_TYPE    ( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_LOWER_MATRIX_TYPE    ( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type5                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type5 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( Type5* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsStrictlyUpper' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsStrictlyUpper' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsStrictlyUpper()
{
   using blaze::DynamicMatrix;
   using blaze::StrictlyLowerMatrix;
   using blaze::UpperMatrix;
   using blaze::StrictlyUpperMatrix;
   using blaze::DiagonalMatrix;

   using Type1 = DynamicMatrix<int>;
   using Type2 = StrictlyLowerMatrix< DynamicMatrix<int> >;
   using Type3 = UpperMatrix< DynamicMatrix<int> >;
   using Type4 = StrictlyUpperMatrix< DynamicMatrix<int> >;
   using Type5 = DiagonalMatrix< DynamicMatrix<int> >;

   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_UPPER_MATRIX_TYPE    ( Type4                 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_UPPER_MATRIX_TYPE    ( Type4 const           );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_UPPER_MATRIX_TYPE    ( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_UPPER_MATRIX_TYPE    ( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type5                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type5 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( Type1* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsSymmetric' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsSymmetric' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsSymmetric()
{
   using blaze::complex;
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::HermitianMatrix;
   using blaze::DiagonalMatrix;
   using blaze::IdentityMatrix;

   using Type1 = DynamicMatrix<int>;
   using Type2 = SymmetricMatrix< DynamicMatrix<int> >;
   using Type3 = HermitianMatrix< DynamicMatrix<int> >;
   using Type4 = HermitianMatrix< DynamicMatrix< complex<int> > >;
   using Type5 = DiagonalMatrix< DynamicMatrix<int> >;
   using Type6 = IdentityMatrix<int>;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type2                 );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type2 const           );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type3                 );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type3 const           );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type5                 );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type5 const           );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type6                 );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type6 const           );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE    ( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( Type6* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsUniform' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsUniform' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsUniform()
{
   using blaze::DynamicVector;
   using blaze::DynamicMatrix;
   using blaze::UniformVector;
   using blaze::UniformMatrix;

   using Type1 = DynamicVector<int>;
   using Type2 = DynamicMatrix<int>;
   using Type3 = UniformVector<int>;
   using Type4 = UniformMatrix<int>;

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_UNIFORM_TYPE    ( Type3                 );
   BLAZE_CONSTRAINT_MUST_BE_UNIFORM_TYPE    ( Type3 const           );
   BLAZE_CONSTRAINT_MUST_BE_UNIFORM_TYPE    ( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_UNIFORM_TYPE    ( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_UNIFORM_TYPE    ( Type4                 );
   BLAZE_CONSTRAINT_MUST_BE_UNIFORM_TYPE    ( Type4 const           );
   BLAZE_CONSTRAINT_MUST_BE_UNIFORM_TYPE    ( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_UNIFORM_TYPE    ( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( Type4* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsUniLower' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsUniLower' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsUniLower()
{
   using blaze::DynamicMatrix;
   using blaze::LowerMatrix;
   using blaze::UniLowerMatrix;
   using blaze::UniUpperMatrix;
   using blaze::DiagonalMatrix;
   using blaze::IdentityMatrix;

   using Type1 = DynamicMatrix<int>;
   using Type2 = LowerMatrix< DynamicMatrix<int> >;
   using Type3 = UniLowerMatrix< DynamicMatrix<int> >;
   using Type4 = UniUpperMatrix< DynamicMatrix<int> >;
   using Type5 = DiagonalMatrix< DynamicMatrix<int> >;
   using Type6 = IdentityMatrix<int>;

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_UNILOWER_MATRIX_TYPE    ( Type3                 );
   BLAZE_CONSTRAINT_MUST_BE_UNILOWER_MATRIX_TYPE    ( Type3 const           );
   BLAZE_CONSTRAINT_MUST_BE_UNILOWER_MATRIX_TYPE    ( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_UNILOWER_MATRIX_TYPE    ( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type5                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type5 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_UNILOWER_MATRIX_TYPE    ( Type6                 );
   BLAZE_CONSTRAINT_MUST_BE_UNILOWER_MATRIX_TYPE    ( Type6 const           );
   BLAZE_CONSTRAINT_MUST_BE_UNILOWER_MATRIX_TYPE    ( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_UNILOWER_MATRIX_TYPE    ( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE( Type6* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsUniUpper' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsUniUpper' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsUniUpper()
{
   using blaze::DynamicMatrix;
   using blaze::UniLowerMatrix;
   using blaze::UniUpperMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;
   using blaze::IdentityMatrix;

   using Type1 = DynamicMatrix<int>;
   using Type2 = UniLowerMatrix< DynamicMatrix<int> >;
   using Type3 = UpperMatrix< DynamicMatrix<int> >;
   using Type4 = UniUpperMatrix< DynamicMatrix<int> >;
   using Type5 = DiagonalMatrix< DynamicMatrix<int> >;
   using Type6 = IdentityMatrix<int>;

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_UNIUPPER_MATRIX_TYPE    ( Type4                 );
   BLAZE_CONSTRAINT_MUST_BE_UNIUPPER_MATRIX_TYPE    ( Type4 const           );
   BLAZE_CONSTRAINT_MUST_BE_UNIUPPER_MATRIX_TYPE    ( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_UNIUPPER_MATRIX_TYPE    ( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type5                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type5 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_UNIUPPER_MATRIX_TYPE    ( Type6                 );
   BLAZE_CONSTRAINT_MUST_BE_UNIUPPER_MATRIX_TYPE    ( Type6 const           );
   BLAZE_CONSTRAINT_MUST_BE_UNIUPPER_MATRIX_TYPE    ( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_UNIUPPER_MATRIX_TYPE    ( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( Type6* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsUpper' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsUpper' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testIsUpper()
{
   using blaze::DynamicMatrix;
   using blaze::SymmetricMatrix;
   using blaze::LowerMatrix;
   using blaze::UniUpperMatrix;
   using blaze::StrictlyUpperMatrix;
   using blaze::UpperMatrix;
   using blaze::DiagonalMatrix;

   using Type1 = DynamicMatrix<int>;
   using Type2 = SymmetricMatrix< DynamicMatrix<int> >;
   using Type3 = LowerMatrix< DynamicMatrix<int> >;
   using Type4 = UpperMatrix< DynamicMatrix<int> >;
   using Type5 = UniUpperMatrix< DynamicMatrix<int> >;
   using Type6 = StrictlyUpperMatrix< DynamicMatrix<int> >;
   using Type7 = DiagonalMatrix< DynamicMatrix<int> >;

   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type4                 );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type4 const           );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type5                 );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type5 const           );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type6                 );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type6 const           );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type6* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type7                 );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type7 const           );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type7 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_UPPER_MATRIX_TYPE    ( Type7 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type7&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type7*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type7* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type7* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE( Type7* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsVector' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsVector' type trait. In case
// an error is detected, a compilation error is created.
*/
void OperationTest::testIsVector()
{
   using blaze::columnVector;
   using blaze::rowVector;
   using blaze::rowMajor;
   using blaze::columnMajor;

   using Type1 = blaze::DynamicVector<float,columnVector>;
   using Type2 = blaze::CompressedVector<int,rowVector>;
   using Type3 = int;
   using Type4 = blaze::complex<float>;
   using Type5 = blaze::DynamicMatrix<double,rowMajor>;
   using Type6 = blaze::CompressedMatrix<int,columnMajor>;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE    ( Type1                 );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE    ( Type1 const           );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE    ( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE    ( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE    ( Type2                 );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE    ( Type2 const           );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE    ( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_TYPE    ( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type3                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type3 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type4                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type4 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type4* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type5                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type5 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type5 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type5 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type5&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type5*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type5* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type5* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type5* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type6                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type6 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type6 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type6 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type6&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type6*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type6* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type6* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VECTOR_TYPE( Type6* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'IsZero' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'IsZero' type trait. In case
// an error is detected, a compilation error is created.
*/
void OperationTest::testIsZero()
{
   using blaze::DynamicVector;
   using blaze::DynamicMatrix;
   using blaze::ZeroVector;
   using blaze::ZeroMatrix;

   using Type1 = blaze::DynamicVector<int>;
   using Type2 = blaze::DynamicMatrix<int>;
   using Type3 = blaze::ZeroVector<int>;
   using Type4 = blaze::ZeroMatrix<int>;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type1                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type1 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type1 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type1 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type1&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type1*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type1* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type1* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type1* const volatile );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type2                 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type2 const           );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type2 volatile        );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type2 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type2&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type2*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type2* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type2* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type2* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE    ( Type3                 );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE    ( Type3 const           );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE    ( Type3 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE    ( Type3 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type3&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type3*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type3* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type3* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type3* const volatile );

   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE    ( Type4                 );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE    ( Type4 const           );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE    ( Type4 volatile        );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE    ( Type4 const volatile  );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type4&                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type4*                );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type4* const          );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type4* volatile       );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( Type4* const volatile );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'MakeComplex' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'MakeComplex' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testMakeComplex()
{
   using blaze::MakeComplex_t;
   using blaze::complex;

   using fcplx = complex<float>;
   using dcplx = complex<double>;
   using lcplx = complex<long double>;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<float      >, fcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<double     >, dcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<long double>, lcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<fcplx      >, fcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<dcplx      >, dcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<lcplx      >, lcplx );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const float      >, const fcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const double     >, const dcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const long double>, const lcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const fcplx      >, const fcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const dcplx      >, const dcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const lcplx      >, const lcplx );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<volatile float      >, volatile fcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<volatile double     >, volatile dcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<volatile long double>, volatile lcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<volatile fcplx      >, volatile fcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<volatile dcplx      >, volatile dcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<volatile lcplx      >, volatile lcplx );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const volatile float      >, const volatile fcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const volatile double     >, const volatile dcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const volatile long double>, const volatile lcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const volatile fcplx      >, const volatile fcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const volatile dcplx      >, const volatile dcplx );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( MakeComplex_t<const volatile lcplx      >, const volatile lcplx );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'RemoveAdaptor' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'RemoveAdaptor' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testRemoveAdaptor()
{
   using blaze::DynamicVector;
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;
   using blaze::SymmetricMatrix;
   using blaze::RemoveAdaptor;

   using Source1 = SymmetricMatrix< DynamicMatrix<int> >;
   using Source2 = const SymmetricMatrix< CompressedMatrix<float> >;
   using Source3 = volatile SymmetricMatrix< DynamicMatrix<double> >;
   using Source4 = int;
   using Source5 = const DynamicVector<int>;
   using Source6 = volatile DynamicMatrix<int>;

   using Result1 = DynamicMatrix<int>;
   using Result2 = const CompressedMatrix<float>;
   using Result3 = volatile DynamicMatrix<double>;
   using Result4 = int;
   using Result5 = const DynamicVector<int>;
   using Result6 = volatile DynamicMatrix<int>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source1>::Type, Result1 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source2>::Type, Result2 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source3>::Type, Result3 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source4>::Type, Result4 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source5>::Type, Result5 );
   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( RemoveAdaptor<Source6>::Type, Result6 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'UnderlyingBuiltin' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'UnderlyingBuiltin' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testUnderlyingBuiltin()
{
   using blaze::complex;
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::UnderlyingBuiltin;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<A                >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<A const          >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<A volatile       >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<A const volatile >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<A&               >::Type, A& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<A*               >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<A* const         >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<A* volatile      >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<A* const volatile>::Type, A* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<B                >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<B const          >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<B volatile       >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<B const volatile >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<B&               >::Type, B& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<B*               >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<B* const         >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<B* volatile      >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<B* const volatile>::Type, B* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<C                >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<C const          >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<C volatile       >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<C const volatile >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<C&               >::Type, C& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<C*               >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<C* const         >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<C* volatile      >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<C* const volatile>::Type, C* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<D                >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<D const          >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<D volatile       >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<D const volatile >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<D&               >::Type, D& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<D*               >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<D* const         >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<D* volatile      >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<D* const volatile>::Type, D* );

   using Type1 = double;                                    // Built-in data type
   using Type2 = complex<float>;                            // Complex data type
   using Type3 = std::vector<double>;                       // Container type
   using Type4 = StaticVector<int,3UL>;                     // Vector with built-in element type
   using Type5 = CompressedVector< DynamicVector<float> >;  // Vector with vector element type

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type1                >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type1 const          >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type1 volatile       >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type1 const volatile >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type1&               >::Type, Type1& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type1*               >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type1* const         >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type1* volatile      >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type1* const volatile>::Type, Type1* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type2                >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type2 const          >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type2 volatile       >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type2 const volatile >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type2&               >::Type, Type2& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type2*               >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type2* const         >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type2* volatile      >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type2* const volatile>::Type, Type2* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type3                >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type3 const          >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type3 volatile       >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type3 const volatile >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type3&               >::Type, Type3& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type3*               >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type3* const         >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type3* volatile      >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type3* const volatile>::Type, Type3* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type4                >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type4 const          >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type4 volatile       >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type4 const volatile >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type4&               >::Type, Type4& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type4*               >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type4* const         >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type4* volatile      >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type4* const volatile>::Type, Type4* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type5                >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type5 const          >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type5 volatile       >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type5 const volatile >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type5&               >::Type, Type5& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type5*               >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type5* const         >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type5* volatile      >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingBuiltin<Type5* const volatile>::Type, Type5* );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'UnderlyingElement' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'UnderlyingElement' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testUnderlyingElement()
{
   using blaze::complex;
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::UnderlyingElement;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<A                >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<A const          >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<A volatile       >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<A const volatile >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<A&               >::Type, A& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<A*               >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<A* const         >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<A* volatile      >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<A* const volatile>::Type, A* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<B                >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<B const          >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<B volatile       >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<B const volatile >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<B&               >::Type, B& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<B*               >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<B* const         >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<B* volatile      >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<B* const volatile>::Type, B* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<C                >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<C const          >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<C volatile       >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<C const volatile >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<C&               >::Type, C& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<C*               >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<C* const         >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<C* volatile      >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<C* const volatile>::Type, C* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<D                >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<D const          >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<D volatile       >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<D const volatile >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<D&               >::Type, D& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<D*               >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<D* const         >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<D* volatile      >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<D* const volatile>::Type, D* );

   using Type1 = double;                                    // Built-in data type
   using Type2 = complex<float>;                            // Complex data type
   using Type3 = std::vector<double>;                       // Container type
   using Type4 = StaticVector<int,3UL>;                     // Vector with built-in element type
   using Type5 = CompressedVector< DynamicVector<float> >;  // Vector with vector element type

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type1                >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type1 const          >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type1 volatile       >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type1 const volatile >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type1&               >::Type, Type1& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type1*               >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type1* const         >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type1* volatile      >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type1* const volatile>::Type, Type1* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type2                >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type2 const          >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type2 volatile       >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type2 const volatile >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type2&               >::Type, Type2& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type2*               >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type2* const         >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type2* volatile      >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type2* const volatile>::Type, Type2* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type3                >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type3 const          >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type3 volatile       >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type3 const volatile >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type3&               >::Type, Type3& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type3*               >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type3* const         >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type3* volatile      >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type3* const volatile>::Type, Type3* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type4                >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type4 const          >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type4 volatile       >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type4 const volatile >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type4&               >::Type, Type4& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type4*               >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type4* const         >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type4* volatile      >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type4* const volatile>::Type, Type4* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type5                >::Type, DynamicVector<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type5 const          >::Type, DynamicVector<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type5 volatile       >::Type, DynamicVector<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type5 const volatile >::Type, DynamicVector<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type5&               >::Type, Type5& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type5*               >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type5* const         >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type5* volatile      >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingElement<Type5* const volatile>::Type, Type5* );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'UnderlyingNumeric' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'UnderlyingNumeric' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testUnderlyingNumeric()
{
   using blaze::complex;
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::UnderlyingNumeric;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<A                >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<A const          >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<A volatile       >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<A const volatile >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<A&               >::Type, A& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<A*               >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<A* const         >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<A* volatile      >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<A* const volatile>::Type, A* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<B                >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<B const          >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<B volatile       >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<B const volatile >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<B&               >::Type, B& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<B*               >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<B* const         >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<B* volatile      >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<B* const volatile>::Type, B* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<C                >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<C const          >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<C volatile       >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<C const volatile >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<C&               >::Type, C& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<C*               >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<C* const         >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<C* volatile      >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<C* const volatile>::Type, C* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<D                >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<D const          >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<D volatile       >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<D const volatile >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<D&               >::Type, D& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<D*               >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<D* const         >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<D* volatile      >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<D* const volatile>::Type, D* );

   using Type1 = double;                                    // Built-in data type
   using Type2 = complex<float>;                            // Complex data type
   using Type3 = std::vector<double>;                       // Container type
   using Type4 = StaticVector<int,3UL>;                     // Vector with built-in element type
   using Type5 = CompressedVector< DynamicVector<float> >;  // Vector with vector element type

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type1                >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type1 const          >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type1 volatile       >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type1 const volatile >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type1&               >::Type, Type1& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type1*               >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type1* const         >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type1* volatile      >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type1* const volatile>::Type, Type1* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type2                >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type2 const          >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type2 volatile       >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type2 const volatile >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type2&               >::Type, Type2& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type2*               >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type2* const         >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type2* volatile      >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type2* const volatile>::Type, Type2* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type3                >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type3 const          >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type3 volatile       >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type3 const volatile >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type3&               >::Type, Type3& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type3*               >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type3* const         >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type3* volatile      >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type3* const volatile>::Type, Type3* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type4                >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type4 const          >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type4 volatile       >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type4 const volatile >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type4&               >::Type, Type4& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type4*               >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type4* const         >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type4* volatile      >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type4* const volatile>::Type, Type4* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type5                >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type5 const          >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type5 volatile       >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type5 const volatile >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type5&               >::Type, Type5& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type5*               >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type5* const         >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type5* volatile      >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingNumeric<Type5* const volatile>::Type, Type5* );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'UnderlyingScalar' type trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the mathematical 'UnderlyingScalar' type trait.
// In case an error is detected, a compilation error is created.
*/
void OperationTest::testUnderlyingScalar()
{
   using blaze::complex;
   using blaze::StaticVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::UnderlyingScalar;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<A                >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<A const          >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<A volatile       >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<A const volatile >::Type, A );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<A&               >::Type, A& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<A*               >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<A* const         >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<A* volatile      >::Type, A* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<A* const volatile>::Type, A* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<B                >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<B const          >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<B volatile       >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<B const volatile >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<B&               >::Type, B& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<B*               >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<B* const         >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<B* volatile      >::Type, B* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<B* const volatile>::Type, B* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<C                >::Type, C );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<C const          >::Type, C );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<C volatile       >::Type, C );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<C const volatile >::Type, C );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<C&               >::Type, C& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<C*               >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<C* const         >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<C* volatile      >::Type, C* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<C* const volatile>::Type, C* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<D                >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<D const          >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<D volatile       >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<D const volatile >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<D&               >::Type, D& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<D*               >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<D* const         >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<D* volatile      >::Type, D* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<D* const volatile>::Type, D* );

   using Type1 = double;                                    // Built-in data type
   using Type2 = complex<float>;                            // Complex data type
   using Type3 = std::vector<double>;                       // Container type
   using Type4 = StaticVector<int,3UL>;                     // Vector with built-in element type
   using Type5 = CompressedVector< DynamicVector<float> >;  // Vector with vector element type

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type1                >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type1 const          >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type1 volatile       >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type1 const volatile >::Type, double );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type1&               >::Type, Type1& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type1*               >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type1* const         >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type1* volatile      >::Type, Type1* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type1* const volatile>::Type, Type1* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type2                >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type2 const          >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type2 volatile       >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type2 const volatile >::Type, complex<float> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type2&               >::Type, Type2& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type2*               >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type2* const         >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type2* volatile      >::Type, Type2* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type2* const volatile>::Type, Type2* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type3                >::Type, Type3 );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type3 const          >::Type, Type3 );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type3 volatile       >::Type, Type3 );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type3 const volatile >::Type, Type3 );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type3&               >::Type, Type3& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type3*               >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type3* const         >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type3* volatile      >::Type, Type3* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type3* const volatile>::Type, Type3* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type4                >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type4 const          >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type4 volatile       >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type4 const volatile >::Type, int );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type4&               >::Type, Type4& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type4*               >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type4* const         >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type4* volatile      >::Type, Type4* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type4* const volatile>::Type, Type4* );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type5                >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type5 const          >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type5 volatile       >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type5 const volatile >::Type, float );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type5&               >::Type, Type5& );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type5*               >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type5* const         >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type5* volatile      >::Type, Type5* );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( UnderlyingScalar<Type5* const volatile>::Type, Type5* );
}
//*************************************************************************************************

} // namespace typetraits

} // namespace mathtest

} // namespace blazetest




//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running mathematical type traits operation test..." << std::endl;

   try
   {
      RUN_TYPETRAITS_OPERATION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during mathematical type traits operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
