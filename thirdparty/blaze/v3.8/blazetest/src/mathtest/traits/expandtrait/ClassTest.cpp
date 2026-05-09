//=================================================================================================
/*!
//  \file src/mathtest/traits/expandtrait/ClassTest.cpp
//  \brief Source file for the ExpandTrait class test
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
#include <blaze/math/CustomVector.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/HybridMatrix.h>
#include <blaze/math/HybridVector.h>
#include <blaze/math/InitializerVector.h>
#include <blaze/math/StaticMatrix.h>
#include <blaze/math/StaticVector.h>
#include <blaze/math/traits/ExpandTrait.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/math/UniformMatrix.h>
#include <blaze/math/UniformVector.h>
#include <blaze/math/ZeroVector.h>
#include <blaze/util/Complex.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveCVRef.h>
#include <blazetest/mathtest/traits/expandtrait/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace traits {

namespace expandtrait {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the ExpandTrait class test.
//
// \exception std::runtime_error Error detected.
*/
ClassTest::ClassTest()
{
   testVectorExpansion();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the 'ExpandTrait' class template for vector expansions.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'ExpandTrait' class template for vector
// expansions. In case an error is detected, a compilation error is created.
*/
void ClassTest::testVectorExpansion()
{
   using namespace blaze;


   // StaticVector
   {
      {
         using VT = StaticVector<int,3UL,columnVector>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = StaticVector<int,3UL,rowVector>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = StaticVector<int,3UL,columnVector>;
         using RT = StaticMatrix<int,3UL,5UL,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = StaticVector<int,3UL,rowVector>;
         using RT = StaticMatrix<int,5UL,3UL,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // HybridVector
   {
      {
         using VT = HybridVector<int,3UL,columnVector>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = HybridVector<int,3UL,rowVector>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = HybridVector<int,3UL,columnVector>;
         using RT = HybridMatrix<int,3UL,5UL,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = HybridVector<int,3UL,rowVector>;
         using RT = HybridMatrix<int,5UL,3UL,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // DynamicVector
   {
      {
         using VT = DynamicVector<int,columnVector>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = DynamicVector<int,rowVector>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = DynamicVector<int,columnVector>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = DynamicVector<int,rowVector>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // CustomVector
   {
      {
         using VT = CustomVector<int,unaligned,unpadded,columnVector>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = CustomVector<int,unaligned,unpadded,rowVector>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = CustomVector<int,unaligned,unpadded,columnVector>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = CustomVector<int,unaligned,unpadded,rowVector>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // UniformVector
   {
      {
         using VT = UniformVector<int,columnVector>;
         using RT = UniformMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = UniformVector<int,rowVector>;
         using RT = UniformMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = UniformVector<int,columnVector>;
         using RT = UniformMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = UniformVector<int,rowVector>;
         using RT = UniformMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // InitializerVector
   {
      {
         using VT = InitializerVector<int,columnVector>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = InitializerVector<int,rowVector>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = InitializerVector<int,columnVector>;
         using RT = DynamicMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = InitializerVector<int,rowVector>;
         using RT = DynamicMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // CompressedVector
   {
      {
         using VT = CompressedVector<int,columnVector>;
         using RT = CompressedMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = CompressedVector<int,rowVector>;
         using RT = CompressedMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = CompressedVector<int,columnVector>;
         using RT = CompressedMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = CompressedVector<int,rowVector>;
         using RT = CompressedMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }

   // ZeroVector
   {
      {
         using VT = ZeroVector<int,columnVector>;
         using RT = ZeroMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = ZeroVector<int,rowVector>;
         using RT = ZeroMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand( std::declval<VT>(), std::declval<size_t>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = ZeroVector<int,columnVector>;
         using RT = ZeroMatrix<int,columnMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
      {
         using VT = ZeroVector<int,rowVector>;
         using RT = ZeroMatrix<int,rowMajor>;
         static_assert( IsSame_v< ExpandTrait_t<VT,5UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( expand<5UL>( std::declval<VT>() ) ) >;
         static_assert( StorageOrder_v<Expr> == StorageOrder_v<RT>, "Non-matching storage order detected" );
      }
   }
}
//*************************************************************************************************

} // namespace expandtrait

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
   std::cout << "   Running ExpandTrait class test..." << std::endl;

   try
   {
      RUN_EXPANDTRAIT_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during ExpandTrait class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
