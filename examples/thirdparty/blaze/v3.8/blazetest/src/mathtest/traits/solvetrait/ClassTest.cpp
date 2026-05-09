//=================================================================================================
/*!
//  \file src/mathtest/traits/solvetrait/ClassTest.cpp
//  \brief Source file for the SolveTrait class test
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
#include <blaze/math/Aliases.h>
#include <blaze/math/CustomMatrix.h>
#include <blaze/math/CustomVector.h>
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
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/traits/SolveTrait.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/math/UniformVector.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/util/Complex.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveCVRef.h>
#include <blazetest/mathtest/traits/solvetrait/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace traits {

namespace solvetrait {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the SolveTrait class test.
//
// \exception std::runtime_error Error detected.
*/
ClassTest::ClassTest()
{
   testSingleSolve();
   testMultiSolve();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the 'SolveTrait' class template for single LSE solver operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'SolveTrait' class template for single LSE
// solver operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testSingleSolve()
{
   using namespace blaze;


   // StaticMatrix/...
   {
      // .../StaticVector
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = UniformVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = UniformVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = UniformVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = UniformVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // HybridMatrix/...
   {
      // .../StaticVector
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = UniformVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = UniformVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = UniformVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = UniformVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // DynamicMatrix/...
   {
      // .../StaticVector
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // CustomMatrix/...
   {
      // .../StaticVector
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // InitializerMatrix/...
   {
      // .../StaticVector
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // SymmetricMatrix<DynamicMatrix>/...
   {
      // .../StaticVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // SymmetricMatrix<DynamicMatrix>/... (complex)
   {
      // .../StaticVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StaticVector<complex<double>,4UL,columnVector>;
            using RT = StaticVector<complex<double>,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StaticVector<complex<double>,4UL,columnVector>;
            using RT = StaticVector<complex<double>,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StaticVector<complex<double>,4UL,rowVector>;
            using RT = StaticVector<complex<double>,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StaticVector<complex<double>,4UL,rowVector>;
            using RT = StaticVector<complex<double>,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HybridVector<complex<double>,6UL,columnVector>;
            using RT = HybridVector<complex<double>,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HybridVector<complex<double>,6UL,columnVector>;
            using RT = HybridVector<complex<double>,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HybridVector<complex<double>,6UL,rowVector>;
            using RT = HybridVector<complex<double>,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HybridVector<complex<double>,6UL,rowVector>;
            using RT = HybridVector<complex<double>,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DynamicVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DynamicVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DynamicVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DynamicVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = CustomVector<complex<double>,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = CustomVector<complex<double>,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = CustomVector<complex<double>,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = CustomVector<complex<double>,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniformVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniformVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniformVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniformVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = InitializerVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = InitializerVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = InitializerVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = InitializerVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // HermitianMatrix<DynamicMatrix>/... (symmetric)
   {
      // .../StaticVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // HermitianMatrix<DynamicMatrix>/... (Hermitian)
   {
      // .../StaticVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StaticVector<complex<double>,4UL,columnVector>;
            using RT = StaticVector<complex<double>,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StaticVector<complex<double>,4UL,columnVector>;
            using RT = StaticVector<complex<double>,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StaticVector<complex<double>,4UL,rowVector>;
            using RT = StaticVector<complex<double>,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StaticVector<complex<double>,4UL,rowVector>;
            using RT = StaticVector<complex<double>,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HybridVector<complex<double>,6UL,columnVector>;
            using RT = HybridVector<complex<double>,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HybridVector<complex<double>,6UL,columnVector>;
            using RT = HybridVector<complex<double>,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HybridVector<complex<double>,6UL,rowVector>;
            using RT = HybridVector<complex<double>,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HybridVector<complex<double>,6UL,rowVector>;
            using RT = HybridVector<complex<double>,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DynamicVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DynamicVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DynamicVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DynamicVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = CustomVector<complex<double>,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = CustomVector<complex<double>,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = CustomVector<complex<double>,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = CustomVector<complex<double>,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniformVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniformVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniformVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniformVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = InitializerVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = InitializerVector<complex<double>,columnVector>;
            using RT = DynamicVector<complex<double>,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = InitializerVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = InitializerVector<complex<double>,rowVector>;
            using RT = DynamicVector<complex<double>,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // LowerMatrix<DynamicMatrix>/...
   {
      // .../StaticVector
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // UniLowerMatrix<DynamicMatrix>/...
   {
      // .../StaticVector
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // UpperMatrix<DynamicMatrix>/...
   {
      // .../StaticVector
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // UniUpperMatrix<DynamicMatrix>/...
   {
      // .../StaticVector
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }

   // DiagonalMatrix<DynamicMatrix>/...
   {
      // .../StaticVector
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,columnVector>;
            using RT = StaticVector<double,4UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticVector<double,4UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,columnVector>;
            using RT = HybridVector<double,6UL,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridVector<double,6UL,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'SolveTrait' class template for multi LSE solver operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'SolveTrait' class template for multi LSE
// solver operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testMultiSolve()
{
   using namespace blaze;


   // StaticMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,3UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,3UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,3UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,3UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = StaticMatrix<complex<double>,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<complex<double>,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<complex<double>,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<complex<double>,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = StaticMatrix<complex<double>,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<complex<double>,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<complex<double>,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<complex<double>,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // HybridMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,4UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,4UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,4UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,4UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,4UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,4UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,4UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,4UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = HybridMatrix<double,6UL,4UL,rowMajor>;
            using RT = HybridMatrix<double,5UL,4UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = HybridMatrix<double,6UL,4UL,columnMajor>;
            using RT = HybridMatrix<double,5UL,4UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = HybridMatrix<double,6UL,4UL,rowMajor>;
            using RT = HybridMatrix<double,5UL,4UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = HybridMatrix<double,6UL,4UL,columnMajor>;
            using RT = HybridMatrix<double,5UL,4UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = HybridMatrix<complex<double>,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<complex<double>,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<complex<double>,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<complex<double>,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = HybridMatrix<complex<double>,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<complex<double>,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<complex<double>,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<complex<double>,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<double,5UL,7UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // DynamicMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = DynamicMatrix<complex<double>,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<complex<double>,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<complex<double>,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<complex<double>,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = DynamicMatrix<complex<double>,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<complex<double>,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<complex<double>,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<complex<double>,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<double,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // CustomMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = CustomMatrix<complex<double>,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<complex<double>,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<complex<double>,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<complex<double>,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = CustomMatrix<complex<double>,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<complex<double>,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<complex<double>,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<complex<double>,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // InitializerMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = InitializerMatrix<complex<double>>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<complex<double>>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = InitializerMatrix<complex<double>>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<complex<double>>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<double>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<double>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // SymmetricMatrix<DynamicMatrix>/... (real)
   {
      // .../StaticMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // SymmetricMatrix<DynamicMatrix>/... (complex)
   {
      // .../StaticMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StaticMatrix<complex<double>,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<complex<double>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StaticMatrix<complex<double>,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<complex<double>,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StaticMatrix<complex<double>,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<complex<double>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StaticMatrix<complex<double>,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<complex<double>,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HybridMatrix<complex<double>,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<complex<double>,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HybridMatrix<complex<double>,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<complex<double>,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HybridMatrix<complex<double>,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<complex<double>,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HybridMatrix<complex<double>,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<complex<double>,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DynamicMatrix<complex<double>,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DynamicMatrix<complex<double>,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DynamicMatrix<complex<double>,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DynamicMatrix<complex<double>,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = CustomMatrix<complex<double>,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = CustomMatrix<complex<double>,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = CustomMatrix<complex<double>,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = CustomMatrix<complex<double>,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniformMatrix<complex<double>,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniformMatrix<complex<double>,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniformMatrix<complex<double>,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniformMatrix<complex<double>,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = InitializerMatrix<complex<double>>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = InitializerMatrix<complex<double>>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // HermitianMatrix<DynamicMatrix>/... (symmetric)
   {
      // .../StaticMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // HermitianMatrix<DynamicMatrix>/... (Hermitian)
   {
      // .../StaticMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StaticMatrix<complex<double>,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<complex<double>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StaticMatrix<complex<double>,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<complex<double>,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StaticMatrix<complex<double>,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<complex<double>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StaticMatrix<complex<double>,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<complex<double>,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HybridMatrix<complex<double>,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<complex<double>,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HybridMatrix<complex<double>,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<complex<double>,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HybridMatrix<complex<double>,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<complex<double>,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HybridMatrix<complex<double>,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<complex<double>,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DynamicMatrix<complex<double>,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DynamicMatrix<complex<double>,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DynamicMatrix<complex<double>,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DynamicMatrix<complex<double>,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = CustomMatrix<complex<double>,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = CustomMatrix<complex<double>,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = CustomMatrix<complex<double>,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = CustomMatrix<complex<double>,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniformMatrix<complex<double>,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniformMatrix<complex<double>,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniformMatrix<complex<double>,rowMajor>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniformMatrix<complex<double>,columnMajor>;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = InitializerMatrix<complex<double>>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = InitializerMatrix<complex<double>>;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // LowerMatrix<DynamicMatrix>/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // UniLowerMatrix<DynamicMatrix>/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // UpperMatrix<DynamicMatrix>/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // UniUpperMatrix<DynamicMatrix>/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // DiagonalMatrix<DynamicMatrix>/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,rowMajor>;
            using RT = StaticMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StaticMatrix<double,5UL,7UL,columnMajor>;
            using RT = StaticMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,rowMajor>;
            using RT = HybridMatrix<double,8UL,6UL,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HybridMatrix<double,8UL,6UL,columnMajor>;
            using RT = HybridMatrix<double,8UL,6UL,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<UniformMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = SymmetricMatrix< UniformMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,rowMajor> >;
            using RT = DynamicMatrix<complex<double>,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<double>,columnMajor> >;
            using RT = DynamicMatrix<complex<double>,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< SolveTrait_t<T1,T2>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( solve( std::declval<T1>(), std::declval<T2>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }
}
//*************************************************************************************************

} // namespace solvetrait

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
   std::cout << "   Running SolveTrait class test..." << std::endl;

   try
   {
      RUN_SOLVETRAIT_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during SolveTrait class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
