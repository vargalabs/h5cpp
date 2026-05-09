//=================================================================================================
/*!
//  \file src/mathtest/traits/rowtrait/ClassTest.cpp
//  \brief Source file for the RowTrait class test
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
#include <blaze/math/CustomMatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/InitializerMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StrictlyLowerMatrix.h>
#include <blaze/math/StrictlyUpperMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/math/UniformMatrix.h>
#include <blaze/math/UniLowerMatrix.h>
#include <blaze/math/UniUpperMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/math/Views.h>
#include <blaze/math/ZeroMatrix.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveCVRef.h>
#include <blazetest/mathtest/traits/rowtrait/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace traits {

namespace rowtrait {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the RowTrait class test.
//
// \exception std::runtime_error Error detected.
*/
ClassTest::ClassTest()
{
   testRowOperation();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the 'RowTrait' class template for row operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'RowTrait' class template for row
// operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testRowOperation()
{
   using namespace blaze;


   // StaticMatrix
   {
      {
         using MT = StaticMatrix<int,3UL,5UL,rowMajor>;
         using RT = StaticVector<int,5UL,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StaticMatrix<int,3UL,5UL,columnMajor>;
         using RT = StaticVector<int,5UL,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StaticMatrix<int,3UL,5UL,rowMajor>;
         using RT = StaticVector<int,5UL,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StaticMatrix<int,3UL,5UL,columnMajor>;
         using RT = StaticVector<int,5UL,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HybridMatrix
   {
      {
         using MT = HybridMatrix<int,3UL,5UL,rowMajor>;
         using RT = HybridVector<int,5UL,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HybridMatrix<int,3UL,5UL,columnMajor>;
         using RT = HybridVector<int,5UL,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HybridMatrix<int,3UL,5UL,rowMajor>;
         using RT = HybridVector<int,5UL,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HybridMatrix<int,3UL,5UL,columnMajor>;
         using RT = HybridVector<int,5UL,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // DynamicMatrix
   {
      {
         using MT = DynamicMatrix<int,rowMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = DynamicMatrix<int,columnMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = DynamicMatrix<int,rowMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = DynamicMatrix<int,columnMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CustomMatrix
   {
      {
         using MT = CustomMatrix<int,unaligned,unpadded,rowMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = CustomMatrix<int,unaligned,unpadded,columnMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = CustomMatrix<int,unaligned,unpadded,rowMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = CustomMatrix<int,unaligned,unpadded,columnMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniformMatrix
   {
      {
         using MT = UniformMatrix<int,rowMajor>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniformMatrix<int,columnMajor>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniformMatrix<int,rowMajor>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniformMatrix<int,columnMajor>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // InitializerMatrix
   {
      {
         using MT = InitializerMatrix<int>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = InitializerMatrix<int>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CompressedMatrix
   {
      {
         using MT = CompressedMatrix<int,rowMajor>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = CompressedMatrix<int,columnMajor>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = CompressedMatrix<int,rowMajor>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = CompressedMatrix<int,columnMajor>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // IdentityMatrix
   {
      {
         using MT = IdentityMatrix<int,rowMajor>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = IdentityMatrix<int,columnMajor>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = IdentityMatrix<int,rowMajor>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = IdentityMatrix<int,columnMajor>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // ZeroMatrix
   {
      {
         using MT = ZeroMatrix<int,rowMajor>;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = ZeroMatrix<int,columnMajor>;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = ZeroMatrix<int,rowMajor>;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = ZeroMatrix<int,columnMajor>;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<DynamicMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<DynamicMatrix> (complex)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<UniformMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< UniformMatrix<int,rowMajor> >;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< UniformMatrix<int,columnMajor> >;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< UniformMatrix<int,rowMajor> >;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< UniformMatrix<int,columnMajor> >;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<ZeroMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,rowMajor> >;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,columnMajor> >;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,rowMajor> >;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,columnMajor> >;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HermitianMatrix<DynamicMatrix> (symmetric)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HermitianMatrix<DynamicMatrix> (Hermitian)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // LowerMatrix<DynamicMatrix>
   {
      {
         using MT = LowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = LowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = LowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = LowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniLowerMatrix<DynamicMatrix>
   {
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // StrictlyLowerMatrix<DynamicMatrix>
   {
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UpperMatrix<DynamicMatrix>
   {
      {
         using MT = UpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniUpperMatrix<DynamicMatrix>
   {
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // StrictlyUpperMatrix<DynamicMatrix>
   {
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // DiagonalMatrix<DynamicMatrix>
   {
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row( std::declval<MT>(), 0UL ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< RowTrait_t<MT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( row<0UL>( std::declval<MT>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }
}
//*************************************************************************************************

} // namespace rowtrait

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
   std::cout << "   Running RowTrait class test..." << std::endl;

   try
   {
      RUN_ROWTRAIT_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during RowTrait class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
