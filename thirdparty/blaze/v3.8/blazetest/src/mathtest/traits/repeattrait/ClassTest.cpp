//=================================================================================================
/*!
//  \file src/mathtest/traits/repeattrait/ClassTest.cpp
//  \brief Source file for the RepeatTrait class test
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
#include <utility>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/CustomMatrix.h>
#include <blaze/math/CustomVector.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/HybridVector.h>
#include <blaze/math/InitializerMatrix.h>
#include <blaze/math/InitializerVector.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/StrictlyLowerMatrix.h>
#include <blaze/math/StrictlyUpperMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/traits/RepeatTrait.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/math/UniformMatrix.h>
#include <blaze/math/UniformVector.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/math/ZeroMatrix.h>
#include <blaze/math/ZeroVector.h>
#include <blaze/util/Complex.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveCVRef.h>
#include <blazetest/mathtest/traits/repeattrait/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace traits {

namespace repeattrait {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the RepeatTrait class test.
//
// \exception std::runtime_error Error detected.
*/
ClassTest::ClassTest()
{
   testVectorRepetition();
   testMatrixRepetition();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the 'RepeatTrait' class template for vector repetition.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'RepeatTrait' class template for vector
// repetition. In case an error is detected, a compilation error is created.
*/
void ClassTest::testVectorRepetition()
{
   using namespace blaze;


   // StaticVector
   {
      {
         using VT = StaticVector<int,3UL,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = StaticVector<int,3UL,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = StaticVector<int,3UL,columnVector>;
         using RT = StaticVector<int,15UL,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = StaticVector<int,3UL,rowVector>;
         using RT = StaticVector<int,15UL,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HybridVector
   {
      {
         using VT = HybridVector<int,3UL,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = HybridVector<int,3UL,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = HybridVector<int,3UL,columnVector>;
         using RT = HybridVector<int,15UL,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = HybridVector<int,3UL,rowVector>;
         using RT = HybridVector<int,15UL,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // DynamicVector
   {
      {
         using VT = DynamicVector<int,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = DynamicVector<int,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = DynamicVector<int,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = DynamicVector<int,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CustomVector
   {
      {
         using VT = CustomVector<int,unaligned,unpadded,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CustomVector<int,unaligned,unpadded,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CustomVector<int,unaligned,unpadded,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CustomVector<int,unaligned,unpadded,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniformVector
   {
      {
         using VT = UniformVector<int,columnVector>;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = UniformVector<int,rowVector>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = UniformVector<int,columnVector>;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = UniformVector<int,rowVector>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // InitializerVector
   {
      {
         using VT = InitializerVector<int,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = InitializerVector<int,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = InitializerVector<int,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = InitializerVector<int,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CompressedVector
   {
      {
         using VT = CompressedVector<int,columnVector>;
         using RT = CompressedVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CompressedVector<int,rowVector>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CompressedVector<int,columnVector>;
         using RT = CompressedVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CompressedVector<int,rowVector>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // ZeroVector
   {
      {
         using VT = ZeroVector<int,columnVector>;
         using RT = ZeroVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = ZeroVector<int,rowVector>;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = ZeroVector<int,columnVector>;
         using RT = ZeroVector<int,columnVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = ZeroVector<int,rowVector>;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< RepeatTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<5UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'RepeatTrait' class template for vector repetition.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'RepeatTrait' class template for vector
// repetition. In case an error is detected, a compilation error is created.
*/
void ClassTest::testMatrixRepetition()
{
   using namespace blaze;


   // StaticMatrix
   {
      {
         using MT = StaticMatrix<int,3UL,5UL,rowMajor>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StaticMatrix<int,3UL,5UL,columnMajor>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StaticMatrix<int,3UL,5UL,rowMajor>;
         using RT = StaticMatrix<int,6UL,15UL,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StaticMatrix<int,3UL,5UL,columnMajor>;
         using RT = StaticMatrix<int,6UL,15UL,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // HybridMatrix
   {
      {
         using MT = HybridMatrix<int,3UL,5UL,rowMajor>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HybridMatrix<int,3UL,5UL,columnMajor>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HybridMatrix<int,3UL,5UL,rowMajor>;
         using RT = HybridMatrix<int,6UL,15UL,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HybridMatrix<int,3UL,5UL,columnMajor>;
         using RT = HybridMatrix<int,6UL,15UL,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // DynamicMatrix
   {
      {
         using MT = DynamicMatrix<int,rowMajor>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = DynamicMatrix<int,columnMajor>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = DynamicMatrix<int,rowMajor>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = DynamicMatrix<int,columnMajor>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // CustomMatrix
   {
      {
         using MT = CustomMatrix<int,unaligned,unpadded,rowMajor>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = CustomMatrix<int,unaligned,unpadded,columnMajor>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = CustomMatrix<int,unaligned,unpadded,rowMajor>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = CustomMatrix<int,unaligned,unpadded,columnMajor>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // UniformMatrix
   {
      {
         using MT = UniformMatrix<int,rowMajor>;
         using RT = UniformMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniformMatrix<int,columnMajor>;
         using RT = UniformMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniformMatrix<int,rowMajor>;
         using RT = UniformMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniformMatrix<int,columnMajor>;
         using RT = UniformMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // InitializerMatrix
   {
      {
         using MT = InitializerMatrix<int>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = InitializerMatrix<int>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // CompressedMatrix
   {
      {
         using MT = CompressedMatrix<int,rowMajor>;
         using RT = CompressedMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = CompressedMatrix<int,columnMajor>;
         using RT = CompressedMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = CompressedMatrix<int,rowMajor>;
         using RT = CompressedMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = CompressedMatrix<int,columnMajor>;
         using RT = CompressedMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // IdentityMatrix
   {
      {
         using MT = IdentityMatrix<int,rowMajor>;
         using RT = CompressedMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = IdentityMatrix<int,columnMajor>;
         using RT = CompressedMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = IdentityMatrix<int,rowMajor>;
         using RT = CompressedMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = IdentityMatrix<int,columnMajor>;
         using RT = CompressedMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // ZeroMatrix
   {
      {
         using MT = ZeroMatrix<int,rowMajor>;
         using RT = ZeroMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = ZeroMatrix<int,columnMajor>;
         using RT = ZeroMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = ZeroMatrix<int,rowMajor>;
         using RT = ZeroMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = ZeroMatrix<int,columnMajor>;
         using RT = ZeroMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // SymmetricMatrix (real)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
         static_assert( IsSame_v< RepeatTrait_t<MT,3UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<3UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
         static_assert( IsSame_v< RepeatTrait_t<MT,3UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<3UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // SymmetricMatrix (complex)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicMatrix<complex<int>,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicMatrix<complex<int>,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicMatrix<complex<int>,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicMatrix<complex<int>,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         static_assert( IsSame_v< RepeatTrait_t<MT,3UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<3UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         static_assert( IsSame_v< RepeatTrait_t<MT,3UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<3UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // HermitianMatrix (symmetric)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
         static_assert( IsSame_v< RepeatTrait_t<MT,3UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<3UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
         static_assert( IsSame_v< RepeatTrait_t<MT,3UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<3UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // HermitianMatrix (Hermitian)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicMatrix<complex<int>,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicMatrix<complex<int>,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicMatrix<complex<int>,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicMatrix<complex<int>,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         static_assert( IsSame_v< RepeatTrait_t<MT,3UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<3UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         static_assert( IsSame_v< RepeatTrait_t<MT,3UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<3UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // LowerMatrix
   {
      {
         using MT = LowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = LowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = LowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = LowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // UniLowerMatrix
   {
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // StrictlyLowerMatrix
   {
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // UpperMatrix
   {
      {
         using MT = UpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // UniUpperMatrix
   {
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // StrictlyUpperMatrix
   {
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // DiagonalMatrix
   {
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat( std::declval<MT>(), std::declval<size_t>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< RepeatTrait_t<MT,2UL,3UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( repeat<2UL,3UL>( std::declval<MT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }
}
//*************************************************************************************************

} // namespace repeattrait

} // namespace traits

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
   std::cout << "   Running RepeatTrait class test..." << std::endl;

   try
   {
      RUN_REPEATTRAIT_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during RepeatTrait class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
