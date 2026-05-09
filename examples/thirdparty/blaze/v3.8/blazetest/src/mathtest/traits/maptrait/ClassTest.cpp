//=================================================================================================
/*!
//  \file src/mathtest/traits/maptrait/ClassTest.cpp
//  \brief Source file for the MapTrait class test
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/HybridVector.h>
#include <blaze/math/IdentityMatrix.h>
#include <blaze/math/InitializerMatrix.h>
#include <blaze/math/InitializerVector.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/StrictlyLowerMatrix.h>
#include <blaze/math/StrictlyUpperMatrix.h>
#include <blaze/math/SymmetricMatrix.h>
#include <blaze/math/traits/MapTrait.h>
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
#include <blazetest/mathtest/traits/maptrait/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace traits {

namespace maptrait {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the MapTrait class test.
//
// \exception std::runtime_error Error detected.
*/
ClassTest::ClassTest()
{
   testUnaryScalarOperation();
   testUnaryVectorOperation();
   testUnaryMatrixOperation();

   testBinaryScalarOperation();
   testBinaryVectorOperation();
   testBinaryMatrixOperation();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the 'MapTrait' class template for unary scalar operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'MapTrait' class template for unary scalar
// operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testUnaryScalarOperation()
{
   using namespace blaze;
   using OP = Pow2;


   // int
   {
      using ST = int;
      using RT = int;
      static_assert( IsSame_v< MapTrait_t<ST,OP>, RT >, "Non-matching type detected" );
   }

   // double
   {
      using ST = double;
      using RT = double;
      static_assert( IsSame_v< MapTrait_t<ST,OP>, RT >, "Non-matching type detected" );
   }

   // complex<int>
   {
      using ST = complex<int>;
      using RT = complex<int>;
      static_assert( IsSame_v< MapTrait_t<ST,OP>, RT >, "Non-matching type detected" );
   }

   // complex<double>
   {
      using ST = complex<double>;
      using RT = complex<double>;
      static_assert( IsSame_v< MapTrait_t<ST,OP>, RT >, "Non-matching type detected" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'MapTrait' class template for unary vector operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'MapTrait' class template for unary vector
// operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testUnaryVectorOperation()
{
   using namespace blaze;
   using OP = Pow2;


   // StaticVector
   {
      {
         using VT = StaticVector<int,3UL,columnVector>;
         using RT = StaticVector<int,3UL,columnVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = StaticVector<int,3UL,rowVector>;
         using RT = StaticVector<int,3UL,rowVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HybridVector
   {
      {
         using VT = HybridVector<int,5UL,columnVector>;
         using RT = HybridVector<int,5UL,columnVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = HybridVector<int,5UL,rowVector>;
         using RT = HybridVector<int,5UL,rowVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // DynamicVector
   {
      {
         using VT = DynamicVector<int,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = DynamicVector<int,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CustomVector
   {
      {
         using VT = CustomVector<int,unaligned,unpadded,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CustomVector<int,unaligned,unpadded,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniformVector
   {
      {
         using VT = UniformVector<int,columnVector>;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = UniformVector<int,rowVector>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // InitializerVector
   {
      {
         using VT = InitializerVector<int,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = InitializerVector<int,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CompressedVector
   {
      {
         using VT = CompressedVector<int,columnVector>;
         using RT = CompressedVector<int,columnVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CompressedVector<int,rowVector>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // ZeroVector
   {
      {
         using VT = ZeroVector<int,columnVector>;
         using RT = ZeroVector<int,columnVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = ZeroVector<int,rowVector>;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< MapTrait_t<VT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<VT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'MapTrait' class template for unary matrix operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'MapTrait' class template for unary matrix
// operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testUnaryMatrixOperation()
{
   using namespace blaze;
   using OP = Pow2;


   // StaticMatrix
   {
      {
         using MT = StaticMatrix<int,3UL,5UL,rowMajor>;
         using RT = StaticMatrix<int,3UL,5UL,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StaticMatrix<int,3UL,5UL,columnMajor>;
         using RT = StaticMatrix<int,3UL,5UL,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // HybridMatrix
   {
      {
         using MT = HybridMatrix<int,5UL,7UL,rowMajor>;
         using RT = HybridMatrix<int,5UL,7UL,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HybridMatrix<int,5UL,7UL,columnMajor>;
         using RT = HybridMatrix<int,5UL,7UL,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // DynamicMatrix
   {
      {
         using MT = DynamicMatrix<int,rowMajor>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = DynamicMatrix<int,columnMajor>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // CustomMatrix
   {
      {
         using MT = CustomMatrix<int,unaligned,unpadded,rowMajor>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = CustomMatrix<int,unaligned,unpadded,columnMajor>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // UniformMatrix
   {
      {
         using MT = UniformMatrix<int,rowMajor>;
         using RT = UniformMatrix<int,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniformMatrix<int,columnMajor>;
         using RT = UniformMatrix<int,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // InitializerMatrix
   {
      {
         using MT = InitializerMatrix<int>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // CompressedMatrix
   {
      {
         using MT = CompressedMatrix<int,rowMajor>;
         using RT = CompressedMatrix<int,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = CompressedMatrix<int,columnMajor>;
         using RT = CompressedMatrix<int,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // IdentityMatrix
   {
      {
         using MT = IdentityMatrix<int,rowMajor>;
         using RT = IdentityMatrix<int,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = IdentityMatrix<int,columnMajor>;
         using RT = IdentityMatrix<int,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // ZeroMatrix
   {
      {
         using MT = ZeroMatrix<int,rowMajor>;
         using RT = ZeroMatrix<int,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = ZeroMatrix<int,columnMajor>;
         using RT = ZeroMatrix<int,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // SymmetricMatrix<DynamicMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // SymmetricMatrix<DynamicMatrix> (complex)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicMatrix<complex<int>,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicMatrix<complex<int>,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // SymmetricMatrix<UniformMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< UniformMatrix<int,rowMajor> >;
         using RT = UniformMatrix<int,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< UniformMatrix<int,columnMajor> >;
         using RT = UniformMatrix<int,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // SymmetricMatrix<ZeroMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,rowMajor> >;
         using RT = ZeroMatrix<int,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,columnMajor> >;
         using RT = ZeroMatrix<int,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // HermitianMatrix<DynamicMatrix> (symmetric)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // HermitianMatrix<DynamicMatrix> (Hermitian)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicMatrix<complex<int>,rowMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicMatrix<complex<int>,columnMajor>;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // LowerMatrix<DynamicMatrix>
   {
      {
         using MT = LowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = LowerMatrix< DynamicMatrix<int,rowMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = LowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = LowerMatrix< DynamicMatrix<int,columnMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // UniLowerMatrix<DynamicMatrix>
   {
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // StrictlyLowerMatrix<DynamicMatrix>
   {
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // UpperMatrix<DynamicMatrix>
   {
      {
         using MT = UpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = UpperMatrix< DynamicMatrix<int,rowMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = UpperMatrix< DynamicMatrix<int,columnMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // UniUpperMatrix<DynamicMatrix>
   {
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // StrictlyUpperMatrix<DynamicMatrix>
   {
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // DiagonalMatrix<DynamicMatrix>
   {
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
         static_assert( IsSame_v< MapTrait_t<MT,OP>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( map( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'MapTrait' class template for binary scalar operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'MapTrait' class template for binary scalar
// operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testBinaryScalarOperation()
{
   using namespace blaze;
   using OP = Mult;


   // int/...
   {
      // .../int
      {
         using T1 = int;
         using T2 = int;
         using RT = int;
         static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );
      }

      // .../double
      {
         using T1 = int;
         using T2 = double;
         using RT = double;
         static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );
      }
   }

   // double
   {
      // .../int
      {
         using T1 = double;
         using T2 = int;
         using RT = double;
         static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );
      }

      // .../double
      {
         using T1 = double;
         using T2 = double;
         using RT = double;
         static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );
      }

      // .../complex<double>
      {
         using T1 = double;
         using T2 = complex<double>;
         using RT = complex<double>;
         static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );
      }
   }

   // complex<double>
   {
      // .../double
      {
         using T1 = complex<double>;
         using T2 = double;
         using RT = complex<double>;
         static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );
      }

      // .../complex<double>
      {
         using T1 = complex<double>;
         using T2 = complex<double>;
         using RT = complex<double>;
         static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'MapTrait' class template for binary vector operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'MapTrait' class template for binary vector
// operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testBinaryVectorOperation()
{
   using namespace blaze;
   using OP = Mult;


   // StaticVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = StaticMatrix<double,3UL,4UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = HybridVector<double,5UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = HybridVector<double,5UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridMatrix<double,3UL,6UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = StaticVector<int,3UL,columnVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // HybridVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = HybridMatrix<double,5UL,4UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = HybridVector<double,5UL,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = HybridVector<double,5UL,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = HybridMatrix<double,5UL,6UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = HybridVector<int,5UL,columnVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // DynamicVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = HybridVector<double,5UL,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = HybridVector<double,5UL,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = DynamicVector<int,columnVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // CustomVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = HybridVector<double,5UL,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = HybridVector<double,5UL,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = CustomVector<int,unaligned,unpadded,columnVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // UniformVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = HybridVector<double,5UL,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = HybridVector<double,5UL,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = UniformVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = UniformVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = UniformMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = UniformVector<int,columnVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // InitializerVector/...
   {
      // .../StaticVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = StaticVector<double,3UL,columnVector>;
            using RT = StaticVector<double,3UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = StaticVector<double,3UL,rowVector>;
            using RT = StaticVector<double,3UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = StaticVector<double,4UL,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = HybridVector<double,5UL,columnVector>;
            using RT = HybridVector<double,5UL,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = HybridVector<double,5UL,rowVector>;
            using RT = HybridVector<double,5UL,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = HybridVector<double,6UL,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = DynamicVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = DynamicVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = CustomVector<double,unaligned,unpadded,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = UniformVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = UniformVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerVector
      {
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = InitializerVector<double,columnVector>;
            using RT = DynamicVector<double,columnVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,rowVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicVector<double,rowVector>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
         }
         {
            using T1 = InitializerVector<int,columnVector>;
            using T2 = InitializerVector<double,rowVector>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'MapTrait' class template for binary matrix operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'MapTrait' class template for binary matrix
// operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testBinaryMatrixOperation()
{
   using namespace blaze;
   using OP = Mult;


   // StaticMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = HybridMatrix<double,4UL,6UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = HybridMatrix<double,4UL,6UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = HybridMatrix<double,4UL,6UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = HybridMatrix<double,4UL,6UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = StaticMatrix<int,3UL,5UL,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,5UL,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,7UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = HybridMatrix<double,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = HybridMatrix<complex<int>,5UL,7UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< HybridMatrix<double,5UL,7UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< HybridMatrix<double,5UL,7UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< HybridMatrix<double,5UL,7UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< HybridMatrix<double,5UL,7UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< HybridMatrix<double,5UL,7UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< HybridMatrix<double,5UL,7UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< HybridMatrix<double,5UL,7UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< HybridMatrix<double,5UL,7UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HybridMatrix<int,5UL,7UL,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< HybridMatrix<double,5UL,7UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,8UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DynamicMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,8UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // UniformMatrix/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,8UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = UniformMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = UniformMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = UniformMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = UniformMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,rowMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniformMatrix<int,columnMajor>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = InitializerMatrix<int>;
            using T2 = StaticMatrix<double,3UL,5UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StaticMatrix<double,3UL,5UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,5UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = InitializerMatrix<int>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = InitializerMatrix<int>;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,8UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,8UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HybridMatrix<int,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<complex<int>,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HybridMatrix<int,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<complex<int>,4UL,8UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HybridMatrix<int,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<complex<int>,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HybridMatrix<int,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<complex<int>,4UL,8UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DynamicMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DynamicMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DynamicMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DynamicMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniformMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniformMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniformMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniformMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = InitializerMatrix<int>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = InitializerMatrix<int>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<double,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<double,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,8UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<double,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<double,4UL,8UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DynamicMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DynamicMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,rowMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StaticMatrix<int,3UL,3UL,columnMajor>;
            using RT = StaticMatrix<complex<int>,3UL,3UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HybridMatrix<int,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<complex<int>,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HybridMatrix<int,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<complex<int>,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HybridMatrix<int,4UL,8UL,rowMajor>;
            using RT = HybridMatrix<complex<int>,4UL,8UL,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HybridMatrix<int,4UL,8UL,columnMajor>;
            using RT = HybridMatrix<complex<int>,4UL,8UL,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DynamicMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DynamicMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DynamicMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DynamicMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = CustomMatrix<int,unaligned,unpadded,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniformMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniformMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniformMatrix<int,rowMajor>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniformMatrix<int,columnMajor>;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = InitializerMatrix<int>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = InitializerMatrix<int>;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DynamicMatrix<complex<int>,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DynamicMatrix<complex<int>,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = LowerMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = LowerMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = LowerMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = LowerMatrix< HybridMatrix<double,4UL,8UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = LowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = LowerMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = LowerMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = LowerMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = LowerMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = LowerMatrix< HybridMatrix<double,4UL,8UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = IdentityMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = IdentityMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = IdentityMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = IdentityMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // StrictlyLowerMatrix<DynamicMatrix>/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StrictlyLowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StrictlyLowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StrictlyLowerMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StrictlyLowerMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = StrictlyLowerMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = StrictlyLowerMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = StrictlyLowerMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = StrictlyLowerMatrix< HybridMatrix<double,4UL,8UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = UpperMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = UpperMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = UpperMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = UpperMatrix< HybridMatrix<double,4UL,8UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = UpperMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = UpperMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = UpperMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = UpperMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = UpperMatrix< HybridMatrix<double,4UL,8UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = IdentityMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = IdentityMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = IdentityMatrix<double,rowMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = IdentityMatrix<double,columnMajor>;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }

   // StrictlyUpperMatrix<DynamicMatrix>/...
   {
      // .../StaticMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StrictlyUpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StrictlyUpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = StrictlyUpperMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = StrictlyUpperMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = StrictlyUpperMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = StrictlyUpperMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = StrictlyUpperMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = StrictlyUpperMatrix< HybridMatrix<double,4UL,8UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
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
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = DiagonalMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = DiagonalMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,rowMajor>;
            using RT = DiagonalMatrix< StaticMatrix<double,3UL,3UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StaticMatrix<double,3UL,3UL,columnMajor>;
            using RT = DiagonalMatrix< StaticMatrix<double,3UL,3UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HybridMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = DiagonalMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = DiagonalMatrix< HybridMatrix<double,4UL,8UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,rowMajor>;
            using RT = DiagonalMatrix< HybridMatrix<double,4UL,8UL,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HybridMatrix<double,4UL,8UL,columnMajor>;
            using RT = DiagonalMatrix< HybridMatrix<double,4UL,8UL,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DynamicMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,rowMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DynamicMatrix<double,columnMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../CustomMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,rowMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = CustomMatrix<double,unaligned,unpadded,columnMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniformMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,rowMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniformMatrix<double,columnMajor>;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../InitializerMatrix
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = InitializerMatrix<double>;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (real)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../SymmetricMatrix<DynamicMatrix> (complex)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (symmetric)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../HermitianMatrix<DynamicMatrix> (Hermitian)
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<complex<int>,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../LowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = LowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyLowerMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyLowerMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../UniUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = UniUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../StrictlyUpperMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = StrictlyUpperMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }

      // .../DiagonalMatrix<DynamicMatrix>
      {
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,rowMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
         {
            using T1 = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
            using T2 = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            using RT = DiagonalMatrix< DynamicMatrix<double,columnMajor> >;
            static_assert( IsSame_v< MapTrait_t<T1,T2,OP>, RT >, "Non-matching type detected" );

            using Expr = RemoveCVRef_t< decltype( map( std::declval<T1>(), std::declval<T2>(), std::declval<OP>() ) ) >;
            static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
            static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
         }
      }
   }
}
//*************************************************************************************************

} // namespace maptrait

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
   std::cout << "   Running MapTrait class test..." << std::endl;

   try
   {
      RUN_MAPTRAIT_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during MapTrait class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
