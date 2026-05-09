//=================================================================================================
/*!
//  \file src/utiltest/valuetraits/OperationTest.cpp
//  \brief Source file for the value traits operation test
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
#include <blaze/util/valuetraits/IsEven.h>
#include <blaze/util/valuetraits/IsMultipleOf.h>
#include <blaze/util/valuetraits/IsOdd.h>
#include <blaze/util/valuetraits/IsPowerOf.h>
#include <blaze/util/StaticAssert.h>
#include <blazetest/utiltest/valuetraits/OperationTest.h>


namespace blazetest {

namespace utiltest {

namespace valuetraits {

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
   testIsEven();
   testIsOdd();
   testIsMultipleOf();
   testIsPowerOf();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST VALUE TRAITS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the 'IsEven' value trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'IsEven' value trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsEven()
{
   using blaze::IsEven;

   BLAZE_STATIC_ASSERT( IsEven<0>::value == 1 );
   BLAZE_STATIC_ASSERT( IsEven<1>::value == 0 );
   BLAZE_STATIC_ASSERT( IsEven<2>::value == 1 );
   BLAZE_STATIC_ASSERT( IsEven<3>::value == 0 );
   BLAZE_STATIC_ASSERT( IsEven<4>::value == 1 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'IsOdd' value trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'IsOdd' value trait. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIsOdd()
{
   using blaze::IsOdd;

   BLAZE_STATIC_ASSERT( IsOdd<0>::value == 0 );
   BLAZE_STATIC_ASSERT( IsOdd<1>::value == 1 );
   BLAZE_STATIC_ASSERT( IsOdd<2>::value == 0 );
   BLAZE_STATIC_ASSERT( IsOdd<3>::value == 1 );
   BLAZE_STATIC_ASSERT( IsOdd<4>::value == 0 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'IsMultipleOf' value trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'IsMultipleOf' value trait. In case an
// error is detected, a compilation error is created.
*/
void OperationTest::testIsMultipleOf()
{
   using blaze::IsMultipleOf;

   BLAZE_STATIC_ASSERT( ( IsMultipleOf<8,2>::value == 1 ) );
   BLAZE_STATIC_ASSERT( ( IsMultipleOf<2,2>::value == 1 ) );
   BLAZE_STATIC_ASSERT( ( IsMultipleOf<0,2>::value == 1 ) );
   BLAZE_STATIC_ASSERT( ( IsMultipleOf<0,0>::value == 1 ) );
   BLAZE_STATIC_ASSERT( ( IsMultipleOf<5,3>::value == 0 ) );
   BLAZE_STATIC_ASSERT( ( IsMultipleOf<2,3>::value == 0 ) );
   BLAZE_STATIC_ASSERT( ( IsMultipleOf<2,0>::value == 0 ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the 'IsPowerOf' value trait.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a compile time test of the 'IsPowerOf' value trait. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testIsPowerOf()
{
   using blaze::IsPowerOf;

   BLAZE_STATIC_ASSERT( ( IsPowerOf<2,8>::value  == 1 ) );
   BLAZE_STATIC_ASSERT( ( IsPowerOf<3,27>::value == 1 ) );
   BLAZE_STATIC_ASSERT( ( IsPowerOf<5,1>::value  == 1 ) );
   BLAZE_STATIC_ASSERT( ( IsPowerOf<1,1>::value  == 1 ) );
   BLAZE_STATIC_ASSERT( ( IsPowerOf<0,0>::value  == 1 ) );
   BLAZE_STATIC_ASSERT( ( IsPowerOf<2,14>::value == 0 ) );
   BLAZE_STATIC_ASSERT( ( IsPowerOf<1,5>::value  == 0 ) );
   BLAZE_STATIC_ASSERT( ( IsPowerOf<0,5>::value  == 0 ) );
   BLAZE_STATIC_ASSERT( ( IsPowerOf<2,0>::value  == 0 ) );
}
//*************************************************************************************************

} // namespace valuetraits

} // namespace utiltest

} // namespace blazetest




//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running value traits operation test..." << std::endl;

   try
   {
      RUN_VALUETRAITS_OPERATION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during value traits operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
