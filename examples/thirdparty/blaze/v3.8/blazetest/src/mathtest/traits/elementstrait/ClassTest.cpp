//=================================================================================================
/*!
//  \file src/mathtest/traits/elementstrait/ClassTest.cpp
//  \brief Source file for the ElementsTrait class test
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
#include <blaze/math/traits/ElementsTrait.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/math/UniformVector.h>
#include <blaze/math/Views.h>
#include <blaze/math/ZeroVector.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveCVRef.h>
#include <blazetest/mathtest/traits/elementstrait/ClassTest.h>


namespace blazetest {

namespace mathtest {

namespace traits {

namespace elementstrait {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the ElementsTrait class test.
//
// \exception std::runtime_error Error detected.
*/
ClassTest::ClassTest()
{
   testElementsOperation();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the 'ElementsTrait' class template for elements operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'ElementsTrait' class template for elements
// operations. In case an error is detected, a compilation error is created.
*/
void ClassTest::testElementsOperation()
{
   using namespace blaze;


   // StaticVector
   {
      {
         using VT = StaticVector<int,3UL,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = StaticVector<int,3UL,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = StaticVector<int,3UL,columnVector>;
         using RT = StaticVector<int,2UL,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = StaticVector<int,3UL,rowVector>;
         using RT = StaticVector<int,2UL,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // HybridVector
   {
      {
         using VT = HybridVector<int,3UL,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = HybridVector<int,3UL,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = HybridVector<int,3UL,columnVector>;
         using RT = StaticVector<int,2UL,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = HybridVector<int,3UL,rowVector>;
         using RT = StaticVector<int,2UL,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // DynamicVector
   {
      {
         using VT = DynamicVector<int,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = DynamicVector<int,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = DynamicVector<int,columnVector>;
         using RT = StaticVector<int,2UL,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = DynamicVector<int,rowVector>;
         using RT = StaticVector<int,2UL,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CustomVector
   {
      {
         using VT = CustomVector<int,unaligned,unpadded,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CustomVector<int,unaligned,unpadded,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CustomVector<int,unaligned,unpadded,columnVector>;
         using RT = StaticVector<int,2UL,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CustomVector<int,unaligned,unpadded,rowVector>;
         using RT = StaticVector<int,2UL,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // UniformVector
   {
      {
         using VT = UniformVector<int,columnVector>;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = UniformVector<int,rowVector>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = UniformVector<int,columnVector>;
         using RT = UniformVector<int,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = UniformVector<int,rowVector>;
         using RT = UniformVector<int,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // InitializerVector
   {
      {
         using VT = InitializerVector<int,columnVector>;
         using RT = DynamicVector<int,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = InitializerVector<int,rowVector>;
         using RT = DynamicVector<int,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = InitializerVector<int,columnVector>;
         using RT = StaticVector<int,2UL,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = InitializerVector<int,rowVector>;
         using RT = StaticVector<int,2UL,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // CompressedVector
   {
      {
         using VT = CompressedVector<int,columnVector>;
         using RT = CompressedVector<int,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CompressedVector<int,rowVector>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CompressedVector<int,columnVector>;
         using RT = CompressedVector<int,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = CompressedVector<int,rowVector>;
         using RT = CompressedVector<int,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }

   // ZeroVector
   {
      {
         using VT = ZeroVector<int,columnVector>;
         using RT = ZeroVector<int,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = ZeroVector<int,rowVector>;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,0UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements( std::declval<VT>(), { 0UL, 2UL } ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = ZeroVector<int,columnVector>;
         using RT = ZeroVector<int,columnVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
      {
         using VT = ZeroVector<int,rowVector>;
         using RT = ZeroVector<int,rowVector>;
         static_assert( IsSame_v< ElementsTrait_t<VT,2UL>, RT >, "Non-matching type detected" );

         using Expr = RemoveCVRef_t< decltype( elements<0UL,2UL>( std::declval<VT>() ) ) >;
         static_assert( TransposeFlag_v<Expr> == TransposeFlag_v<RT>, "Non-matching transpose flag detected" );
      }
   }
}
//*************************************************************************************************

} // namespace elementstrait

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
   std::cout << "   Running ElementsTrait class test..." << std::endl;

   try
   {
      RUN_ELEMENTSTRAIT_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during ElementsTrait class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
