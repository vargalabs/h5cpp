//=================================================================================================
/*!
//  \file src/mathtest/traits/crosstrait/ClassTest.cpp
//  \brief Source file for the CrossTrait class test
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
#include <blaze/math/CompressedVector.h>
#include <blaze/math/CustomVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/HybridVector.h>
#include <blaze/math/InitializerVector.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/math/UniformVector.h>
#include <blaze/math/ZeroVector.h>
#include <blaze/util/Complex.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveCVRef.h>
#include <blazetest/mathtest/traits/crosstrait/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace traits {

namespace crosstrait {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the CrossTrait class test.
//
// \exception std::runtime_error Error detected.
*/
ClassTest::ClassTest()
{
   testCrossProduct();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the 'CrossTrait' class template for cross product operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'CrossTrait' class template for cross
// product operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testCrossProduct()
{
   using namespace blaze;


   // StaticVector/...
   {
      // .../StaticVector
      {
         using T1 = StaticVector<int,3UL,columnVector>;
         using T2 = StaticVector<double,3UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../HybridVector
      {
         using T1 = StaticVector<int,3UL,columnVector>;
         using T2 = HybridVector<double,7UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../DynamicVector
      {
         using T1 = StaticVector<int,3UL,columnVector>;
         using T2 = DynamicVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CustomVector
      {
         using T1 = StaticVector<int,3UL,columnVector>;
         using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../UniformVector
      {
         using T1 = StaticVector<int,3UL,columnVector>;
         using T2 = UniformVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../InitializerVector
      {
         using T1 = StaticVector<int,3UL,columnVector>;
         using T2 = InitializerVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CompressedVector
      {
         using T1 = StaticVector<int,3UL,columnVector>;
         using T2 = CompressedVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../ZeroVector
      {
         using T1 = StaticVector<int,3UL,columnVector>;
         using T2 = ZeroVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HybridVector/...
   {
      // .../StaticVector
      {
         using T1 = HybridVector<int,5UL,columnVector>;
         using T2 = StaticVector<double,3UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../HybridVector
      {
         using T1 = HybridVector<int,5UL,columnVector>;
         using T2 = HybridVector<double,7UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../DynamicVector
      {
         using T1 = HybridVector<int,5UL,columnVector>;
         using T2 = DynamicVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CustomVector
      {
         using T1 = HybridVector<int,5UL,columnVector>;
         using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../UniformVector
      {
         using T1 = HybridVector<int,5UL,columnVector>;
         using T2 = UniformVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../InitializerVector
      {
         using T1 = HybridVector<int,5UL,columnVector>;
         using T2 = InitializerVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CompressedVector
      {
         using T1 = HybridVector<int,5UL,columnVector>;
         using T2 = CompressedVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../ZeroVector
      {
         using T1 = HybridVector<int,5UL,columnVector>;
         using T2 = ZeroVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // DynamicVector/...
   {
      // .../StaticVector
      {
         using T1 = DynamicVector<int,columnVector>;
         using T2 = StaticVector<double,3UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../HybridVector
      {
         using T1 = DynamicVector<int,columnVector>;
         using T2 = HybridVector<double,7UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../DynamicVector
      {
         using T1 = DynamicVector<int,columnVector>;
         using T2 = DynamicVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CustomVector
      {
         using T1 = DynamicVector<int,columnVector>;
         using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../UniformVector
      {
         using T1 = DynamicVector<int,columnVector>;
         using T2 = UniformVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../InitializerVector
      {
         using T1 = DynamicVector<int,columnVector>;
         using T2 = InitializerVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CompressedVector
      {
         using T1 = DynamicVector<int,columnVector>;
         using T2 = CompressedVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../ZeroVector
      {
         using T1 = DynamicVector<int,columnVector>;
         using T2 = ZeroVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CustomVector/...
   {
      // .../StaticVector
      {
         using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
         using T2 = StaticVector<double,3UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../HybridVector
      {
         using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
         using T2 = HybridVector<double,7UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../DynamicVector
      {
         using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
         using T2 = DynamicVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CustomVector
      {
         using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
         using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../UniformVector
      {
         using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
         using T2 = UniformVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../InitializerVector
      {
         using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
         using T2 = InitializerVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CompressedVector
      {
         using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
         using T2 = CompressedVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../ZeroVector
      {
         using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
         using T2 = ZeroVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniformVector/...
   {
      // .../StaticVector
      {
         using T1 = UniformVector<int,columnVector>;
         using T2 = StaticVector<double,3UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../HybridVector
      {
         using T1 = UniformVector<int,columnVector>;
         using T2 = HybridVector<double,7UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../DynamicVector
      {
         using T1 = UniformVector<int,columnVector>;
         using T2 = DynamicVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CustomVector
      {
         using T1 = UniformVector<int,columnVector>;
         using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../UniformVector
      {
         using T1 = UniformVector<int,columnVector>;
         using T2 = UniformVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../InitializerVector
      {
         using T1 = UniformVector<int,columnVector>;
         using T2 = InitializerVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CompressedVector
      {
         using T1 = UniformVector<int,columnVector>;
         using T2 = CompressedVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../ZeroVector
      {
         using T1 = UniformVector<int,columnVector>;
         using T2 = ZeroVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // InitializerVector/...
   {
      // .../StaticVector
      {
         using T1 = InitializerVector<int,columnVector>;
         using T2 = StaticVector<double,3UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../HybridVector
      {
         using T1 = InitializerVector<int,columnVector>;
         using T2 = HybridVector<double,7UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../DynamicVector
      {
         using T1 = InitializerVector<int,columnVector>;
         using T2 = DynamicVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CustomVector
      {
         using T1 = InitializerVector<int,columnVector>;
         using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../UniformVector
      {
         using T1 = InitializerVector<int,columnVector>;
         using T2 = UniformVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../InitializerVector
      {
         using T1 = InitializerVector<int,columnVector>;
         using T2 = InitializerVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CompressedVector
      {
         using T1 = InitializerVector<int,columnVector>;
         using T2 = CompressedVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../ZeroVector
      {
         using T1 = InitializerVector<int,columnVector>;
         using T2 = ZeroVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CompressedVector/...
   {
      // .../StaticVector
      {
         using T1 = CompressedVector<int,columnVector>;
         using T2 = StaticVector<double,3UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../HybridVector
      {
         using T1 = CompressedVector<int,columnVector>;
         using T2 = HybridVector<double,7UL,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../DynamicVector
      {
         using T1 = CompressedVector<int,columnVector>;
         using T2 = DynamicVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CustomVector
      {
         using T1 = CompressedVector<int,columnVector>;
         using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../UniformVector
      {
         using T1 = CompressedVector<int,columnVector>;
         using T2 = UniformVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../InitializerVector
      {
         using T1 = CompressedVector<int,columnVector>;
         using T2 = InitializerVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CompressedVector
      {
         using T1 = CompressedVector<int,columnVector>;
         using T2 = CompressedVector<double,columnVector>;
         using RT = StaticVector<double,3UL,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../ZeroVector
      {
         using T1 = CompressedVector<int,columnVector>;
         using T2 = ZeroVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // ZeroVector/...
   {
      // .../StaticVector
      {
         using T1 = ZeroVector<int,columnVector>;
         using T2 = StaticVector<double,3UL,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../HybridVector
      {
         using T1 = ZeroVector<int,columnVector>;
         using T2 = HybridVector<double,7UL,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../DynamicVector
      {
         using T1 = ZeroVector<int,columnVector>;
         using T2 = DynamicVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CustomVector
      {
         using T1 = ZeroVector<int,columnVector>;
         using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../UniformVector
      {
         using T1 = ZeroVector<int,columnVector>;
         using T2 = UniformVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../InitializerVector
      {
         using T1 = ZeroVector<int,columnVector>;
         using T2 = InitializerVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../CompressedVector
      {
         using T1 = ZeroVector<int,columnVector>;
         using T2 = CompressedVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }

      // .../ZeroVector
      {
         using T1 = ZeroVector<int,columnVector>;
         using T2 = ZeroVector<double,columnVector>;
         using RT = ZeroVector<double,columnVector>;
         static_assert( IsSame_v< CrossTrait_t<T1,T2>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( std::declval<T1>() % std::declval<T2>() ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }
}
//*************************************************************************************************

} // namespace crosstrait

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
   std::cout << "   Running CrossTrait class test..." << std::endl;

   try
   {
      RUN_CROSSTRAIT_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during CrossTrait class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
