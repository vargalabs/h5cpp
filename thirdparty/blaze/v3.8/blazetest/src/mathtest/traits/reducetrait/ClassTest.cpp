//=================================================================================================
/*!
//  \file src/mathtest/traits/reducetrait/ClassTest.cpp
//  \brief Source file for the ReduceTrait class test
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
#include <blaze/math/traits/ReduceTrait.h>
#include <blaze/math/typetraits/StorageOrder.h>
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
#include <blazetest/mathtest/traits/reducetrait/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace traits {

namespace reducetrait {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the ReduceTrait class test.
//
// \exception std::runtime_error Error detected.
*/
ClassTest::ClassTest()
{
   testTotalVectorReduction();
   testTotalMatrixReduction();
   testRowwiseMatrixReduction();
   testColumnwiseMatrixReduction();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the 'ReduceTrait' class template for total vector reductions.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'ReduceTrait' class template for total
// vector reductions. In case an error is detected, a compilation error is created.
*/
void ClassTest::testTotalVectorReduction()
{
   using namespace blaze;
   using OP = Add;


   // StaticVector
   {
      {
         using VT = StaticVector<int,3UL,columnVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using VT = StaticVector<int,3UL,rowVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // HybridVector
   {
      {
         using VT = HybridVector<int,3UL,columnVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using VT = HybridVector<int,3UL,rowVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // DynamicVector
   {
      {
         using VT = DynamicVector<int,columnVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using VT = DynamicVector<int,rowVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // CustomVector
   {
      {
         using VT = CustomVector<int,unaligned,unpadded,columnVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using VT = CustomVector<int,unaligned,unpadded,rowVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // UniformVector
   {
      {
         using VT = UniformVector<int,columnVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using VT = UniformVector<int,rowVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // InitializerVector
   {
      {
         using VT = InitializerVector<int,columnVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using VT = InitializerVector<int,rowVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // CompressedVector
   {
      {
         using VT = CompressedVector<int,columnVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using VT = CompressedVector<int,rowVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // ZeroVector
   {
      {
         using VT = ZeroVector<int,columnVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using VT = ZeroVector<int,rowVector>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<VT,OP>, RT >, "Non-matching type detected" );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'ReduceTrait' class template for total matrix reductions.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'ReduceTrait' class template for total
// matrix reductions. In case an error is detected, a compilation error is created.
*/
void ClassTest::testTotalMatrixReduction()
{
   using namespace blaze;
   using OP = Add;


   // StaticMatrix
   {
      {
         using MT = StaticMatrix<int,3UL,5UL,rowMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = StaticMatrix<int,3UL,5UL,columnMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // HybridMatrix
   {
      {
         using MT = HybridMatrix<int,3UL,5UL,rowMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = HybridMatrix<int,3UL,5UL,columnMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // DynamicMatrix
   {
      {
         using MT = DynamicMatrix<int,rowMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = DynamicMatrix<int,columnMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // CustomMatrix
   {
      {
         using MT = CustomMatrix<int,unaligned,unpadded,rowMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = CustomMatrix<int,unaligned,unpadded,columnMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // UniformMatrix
   {
      {
         using MT = UniformMatrix<int,rowMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = UniformMatrix<int,columnMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // InitializerMatrix
   {
      {
         using MT = InitializerMatrix<int>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // CompressedMatrix
   {
      {
         using MT = CompressedMatrix<int,rowMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = CompressedMatrix<int,columnMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // IdentityMatrix
   {
      {
         using MT = IdentityMatrix<int,rowMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = IdentityMatrix<int,columnMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // ZeroMatrix
   {
      {
         using MT = ZeroMatrix<int,rowMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = ZeroMatrix<int,columnMajor>;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // SymmetricMatrix<DynamicMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // SymmetricMatrix<DynamicMatrix> (complex)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = complex<int>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = complex<int>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // SymmetricMatrix<UniformMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< UniformMatrix<int,rowMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = SymmetricMatrix< UniformMatrix<int,columnMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // SymmetricMatrix<ZeroMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,rowMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,columnMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // HermitianMatrix<DynamicMatrix> (symmetric)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // HermitianMatrix<DynamicMatrix> (Hermitian)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = complex<int>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = complex<int>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // LowerMatrix<DynamicMatrix>
   {
      {
         using MT = LowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = LowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // UniLowerMatrix<DynamicMatrix>
   {
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // StrictlyLowerMatrix<DynamicMatrix>
   {
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // UpperMatrix<DynamicMatrix>
   {
      {
         using MT = UpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = UpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // UniUpperMatrix<DynamicMatrix>
   {
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // StrictlyUpperMatrix<DynamicMatrix>
   {
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }

   // DiagonalMatrix<DynamicMatrix>
   {
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = int;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP>, RT >, "Non-matching type detected" );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'ReduceTrait' class template for rowwise matrix reductions.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'ReduceTrait' class template for rowwise
// matrix reductions. In case an error is detected, a compilation error is created.
*/
void ClassTest::testRowwiseMatrixReduction()
{
   using namespace blaze;
   using OP = Add;


   // StaticMatrix
   {
      {
         using MT = StaticMatrix<int,3UL,5UL,rowMajor>;
         using RT = StaticVector<int,3UL,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StaticMatrix<int,3UL,5UL,columnMajor>;
         using RT = StaticVector<int,3UL,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HybridMatrix
   {
      {
         using MT = HybridMatrix<int,3UL,5UL,rowMajor>;
         using RT = HybridVector<int,3UL,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HybridMatrix<int,3UL,5UL,columnMajor>;
         using RT = HybridVector<int,3UL,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // DynamicMatrix
   {
      {
         using MT = DynamicMatrix<int,rowMajor>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = DynamicMatrix<int,columnMajor>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CustomMatrix
   {
      {
         using MT = CustomMatrix<int,unaligned,unpadded,rowMajor>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = CustomMatrix<int,unaligned,unpadded,columnMajor>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniformMatrix
   {
      {
         using MT = UniformMatrix<int,rowMajor>;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniformMatrix<int,columnMajor>;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // InitializerMatrix
   {
      {
         using MT = InitializerMatrix<int>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CompressedMatrix
   {
      {
         using MT = CompressedMatrix<int,rowMajor>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = CompressedMatrix<int,columnMajor>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // IdentityMatrix
   {
      {
         using MT = IdentityMatrix<int,rowMajor>;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = IdentityMatrix<int,columnMajor>;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // ZeroMatrix
   {
      {
         using MT = ZeroMatrix<int,rowMajor>;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = ZeroMatrix<int,columnMajor>;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<DynamicMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<DynamicMatrix> (complex)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicVector<complex<int>,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicVector<complex<int>,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<UniformMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< UniformMatrix<int,rowMajor> >;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< UniformMatrix<int,columnMajor> >;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<ZeroMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,rowMajor> >;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,columnMajor> >;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HermitianMatrix<DynamicMatrix> (symmetric)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HermitianMatrix<DynamicMatrix> (Hermitian)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicVector<complex<int>,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicVector<complex<int>,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // LowerMatrix<DynamicMatrix>
   {
      {
         using MT = LowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = LowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniLowerMatrix<DynamicMatrix>
   {
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // StrictlyLowerMatrix<DynamicMatrix>
   {
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UpperMatrix<DynamicMatrix>
   {
      {
         using MT = UpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniUpperMatrix<DynamicMatrix>
   {
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // StrictlyUpperMatrix<DynamicMatrix>
   {
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // DiagonalMatrix<DynamicMatrix>
   {
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,rowwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<rowwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'ReduceTrait' class template for columnwise matrix reductions.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'ReduceTrait' class template for columnwise
// matrix reductions. In case an error is detected, a compilation error is created.
*/
void ClassTest::testColumnwiseMatrixReduction()
{
   using namespace blaze;
   using OP = Add;


   // StaticMatrix
   {
      {
         using MT = StaticMatrix<int,3UL,5UL,rowMajor>;
         using RT = StaticVector<int,5UL,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StaticMatrix<int,3UL,5UL,columnMajor>;
         using RT = StaticVector<int,5UL,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HybridMatrix
   {
      {
         using MT = HybridMatrix<int,3UL,5UL,rowMajor>;
         using RT = HybridVector<int,5UL,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HybridMatrix<int,3UL,5UL,columnMajor>;
         using RT = HybridVector<int,5UL,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // DynamicMatrix
   {
      {
         using MT = DynamicMatrix<int,rowMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = DynamicMatrix<int,columnMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CustomMatrix
   {
      {
         using MT = CustomMatrix<int,unaligned,unpadded,rowMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = CustomMatrix<int,unaligned,unpadded,columnMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniformMatrix
   {
      {
         using MT = UniformMatrix<int,rowMajor>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniformMatrix<int,columnMajor>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // InitializerMatrix
   {
      {
         using MT = DynamicMatrix<int>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CompressedMatrix
   {
      {
         using MT = CompressedMatrix<int,rowMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = CompressedMatrix<int,columnMajor>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // IdentityMatrix
   {
      {
         using MT = IdentityMatrix<int,rowMajor>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = IdentityMatrix<int,columnMajor>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // ZeroMatrix
   {
      {
         using MT = ZeroMatrix<int,rowMajor>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = ZeroMatrix<int,columnMajor>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<DynamicMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<DynamicMatrix> (complex)
   {
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<UniformMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< UniformMatrix<int,rowMajor> >;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< UniformMatrix<int,columnMajor> >;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // SymmetricMatrix<ZeroMatrix> (real)
   {
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,rowMajor> >;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = SymmetricMatrix< ZeroMatrix<int,columnMajor> >;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HermitianMatrix<DynamicMatrix> (symmetric)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HermitianMatrix<DynamicMatrix> (Hermitian)
   {
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,rowMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = HermitianMatrix< DynamicMatrix<complex<int>,columnMajor> >;
         using RT = DynamicVector<complex<int>,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // LowerMatrix<DynamicMatrix>
   {
      {
         using MT = LowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = LowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniLowerMatrix<DynamicMatrix>
   {
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // StrictlyLowerMatrix<DynamicMatrix>
   {
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StrictlyLowerMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UpperMatrix<DynamicMatrix>
   {
      {
         using MT = UpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniUpperMatrix<DynamicMatrix>
   {
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = UniUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // StrictlyUpperMatrix<DynamicMatrix>
   {
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = StrictlyUpperMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // DiagonalMatrix<DynamicMatrix>
   {
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,rowMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using MT = DiagonalMatrix< DynamicMatrix<int,columnMajor> >;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ReduceTrait_t<MT,OP,columnwise>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( reduce<columnwise>( std::declval<MT>(), std::declval<OP>() ) ) >;
         static_assert( IsSame_v< ResultType_t<Expr>, RT >, "Non-matching type detected" );
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }
}
//*************************************************************************************************

} // namespace reducetrait

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
   std::cout << "   Running ReduceTrait class test..." << std::endl;

   try
   {
      RUN_REDUCETRAIT_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during ReduceTrait class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
