//=================================================================================================
/*!
//  \file src/mathtest/shims/OperationTest.cpp
//  \brief Source file for the mathematical shims operation test
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
#include <sstream>
#include <stdexcept>
#include <blaze/math/shims/NextMultiple.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Sign.h>
#include <blazetest/mathtest/shims/OperationTest.h>


namespace blazetest {

namespace mathtest {

namespace shims {

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
   testSign();
   testNextMultiple();
   testPrevMultiple();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the mathematical 'sign' shim.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the mathematical 'sign' shim. In case an error is detected,
// a \a std::runtime_error exception is thrown.
*/
void OperationTest::testSign()
{
   using blaze::sign;

   test_ = "sign() shim";

   // 'int'
   if( sign( 1 ) != 1 || sign( 0 ) != 0 || sign( -1 ) != -1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: sign() with 'int' failed\n"
          << " Details:\n"
          << "   sign(  1 ) = " << sign(  1 ) << " (expected  1)\n"
          << "   sign(  0 ) = " << sign(  0 ) << " (expected  0)\n"
          << "   sign( -1 ) = " << sign( -1 ) << " (expected -1)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'unsigned int'
   if( sign( 1U ) != 1 || sign( 0U ) != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: sign() with 'unsigned int' failed\n"
          << " Details:\n"
          << "   sign( 1U ) = " << sign( 1U ) << " (expected 1)\n"
          << "   sign( 0U ) = " << sign( 0U ) << " (expected 0)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'long'
   if( sign( 1L ) != 1 || sign( 0L ) != 0 || sign( -1L ) != -1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: sign() with 'long' failed\n"
          << " Details:\n"
          << "   sign(  1L ) = " << sign(  1 ) << " (expected  1)\n"
          << "   sign(  0L ) = " << sign(  0 ) << " (expected  0)\n"
          << "   sign( -1L ) = " << sign( -1 ) << " (expected -1)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'unsigned long'
   if( sign( 1UL ) != 1 || sign( 0UL ) != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: sign() with 'unsigned long' failed\n"
          << " Details:\n"
          << "   sign( 1UL ) = " << sign( 1UL ) << " (expected 1)\n"
          << "   sign( 0UL ) = " << sign( 0UL ) << " (expected 0)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'float'
   if( sign( 1.0F ) != 1 || sign( 0.0F ) != 0 || sign( -1.0F ) != -1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: sign() with 'long' failed\n"
          << " Details:\n"
          << "   sign(  1.0F ) = " << sign(  1.0F ) << " (expected  1)\n"
          << "   sign(  0.0F ) = " << sign(  0.0F ) << " (expected  0)\n"
          << "   sign( -1.0F ) = " << sign( -1.0F ) << " (expected -1)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'double'
   if( sign( 1.0 ) != 1 || sign( 0.0 ) != 0 || sign( -1.0 ) != -1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: sign() with 'long' failed\n"
          << " Details:\n"
          << "   sign(  1.0 ) = " << sign(  1.0 ) << " (expected  1)\n"
          << "   sign(  0.0 ) = " << sign(  0.0 ) << " (expected  0)\n"
          << "   sign( -1.0 ) = " << sign( -1.0 ) << " (expected -1)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'long double'
   if( sign( 1.0L ) != 1 || sign( 0.0L ) != 0 || sign( -1.0L ) != -1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: sign() with 'long' failed\n"
          << " Details:\n"
          << "   sign(  1.0L ) = " << sign(  1.0L ) << " (expected  1)\n"
          << "   sign(  0.0L ) = " << sign(  0.0L ) << " (expected  0)\n"
          << "   sign( -1.0L ) = " << sign( -1.0L ) << " (expected -1)\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'nextMultiple' shim.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the mathematical 'nextMultiple' shim. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
void OperationTest::testNextMultiple()
{
   using blaze::nextMultiple;

   test_ = "nextMultiple() shim";

   // 'int'
   if( nextMultiple( 7, 4 ) !=  8 ||
       nextMultiple( 8, 4 ) !=  8 ||
       nextMultiple( 9, 4 ) != 12 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: nextMultiple() with 'int' failed\n"
          << " Details:\n"
          << "   nextMultiple( 7, 4 ) = " << nextMultiple( 7, 4 ) << " (expected  8)\n"
          << "   nextMultiple( 8, 4 ) = " << nextMultiple( 8, 4 ) << " (expected  8)\n"
          << "   nextMultiple( 9, 4 ) = " << nextMultiple( 9, 4 ) << " (expected 12)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'unsigned int'
   if( nextMultiple( 7U, 4U ) !=  8U ||
       nextMultiple( 8U, 4U ) !=  8U ||
       nextMultiple( 9U, 4U ) != 12U ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: nextMultiple() with 'unsigned int' failed\n"
          << " Details:\n"
          << "   nextMultiple( 7, 4 ) = " << nextMultiple( 7U, 4U ) << " (expected  8)\n"
          << "   nextMultiple( 8, 4 ) = " << nextMultiple( 8U, 4U ) << " (expected  8)\n"
          << "   nextMultiple( 9, 4 ) = " << nextMultiple( 9U, 4U ) << " (expected 12)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'long'
   if( nextMultiple( 7L, 4L ) !=  8L ||
       nextMultiple( 8L, 4L ) !=  8L ||
       nextMultiple( 9L, 4L ) != 12L ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: nextMultiple() with 'int' failed\n"
          << " Details:\n"
          << "   nextMultiple( 7, 4 ) = " << nextMultiple( 7L, 4L ) << " (expected  8)\n"
          << "   nextMultiple( 8, 4 ) = " << nextMultiple( 8L, 4L ) << " (expected  8)\n"
          << "   nextMultiple( 9, 4 ) = " << nextMultiple( 9L, 4L ) << " (expected 12)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'unsigned long'
   if( nextMultiple( 7UL, 4UL ) !=  8UL ||
       nextMultiple( 8UL, 4UL ) !=  8UL ||
       nextMultiple( 9UL, 4UL ) != 12UL ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: nextMultiple() with 'unsigned int' failed\n"
          << " Details:\n"
          << "   nextMultiple( 7, 4 ) = " << nextMultiple( 7UL, 4UL ) << " (expected  8)\n"
          << "   nextMultiple( 8, 4 ) = " << nextMultiple( 8UL, 4UL ) << " (expected  8)\n"
          << "   nextMultiple( 9, 4 ) = " << nextMultiple( 9UL, 4UL ) << " (expected 12)\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the mathematical 'prevMultiple' shim.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the mathematical 'prevMultiple' shim. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
void OperationTest::testPrevMultiple()
{
   using blaze::prevMultiple;

   test_ = "prevMultiple() shim";

   // 'int'
   if( prevMultiple( 7, 4 ) != 4 ||
       prevMultiple( 8, 4 ) != 8 ||
       prevMultiple( 9, 4 ) != 8 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: prevMultiple() with 'int' failed\n"
          << " Details:\n"
          << "   prevMultiple( 7, 4 ) = " << prevMultiple( 7, 4 ) << " (expected 4)\n"
          << "   prevMultiple( 8, 4 ) = " << prevMultiple( 8, 4 ) << " (expected 8)\n"
          << "   prevMultiple( 9, 4 ) = " << prevMultiple( 9, 4 ) << " (expected 8)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'unsigned int'
   if( prevMultiple( 7U, 4U ) != 4U ||
       prevMultiple( 8U, 4U ) != 8U ||
       prevMultiple( 9U, 4U ) != 8U ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: prevMultiple() with 'unsigned int' failed\n"
          << " Details:\n"
          << "   prevMultiple( 7, 4 ) = " << prevMultiple( 7U, 4U ) << " (expected 4)\n"
          << "   prevMultiple( 8, 4 ) = " << prevMultiple( 8U, 4U ) << " (expected 8)\n"
          << "   prevMultiple( 9, 4 ) = " << prevMultiple( 9U, 4U ) << " (expected 8)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'long'
   if( prevMultiple( 7L, 4L ) != 4L ||
       prevMultiple( 8L, 4L ) != 8L ||
       prevMultiple( 9L, 4L ) != 8L ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: prevMultiple() with 'int' failed\n"
          << " Details:\n"
          << "   prevMultiple( 7, 4 ) = " << prevMultiple( 7L, 4L ) << " (expected 4)\n"
          << "   prevMultiple( 8, 4 ) = " << prevMultiple( 8L, 4L ) << " (expected 8)\n"
          << "   prevMultiple( 9, 4 ) = " << prevMultiple( 9L, 4L ) << " (expected 8)\n";
      throw std::runtime_error( oss.str() );
   }

   // 'unsigned long'
   if( prevMultiple( 7UL, 4UL ) != 4UL ||
       prevMultiple( 8UL, 4UL ) != 8UL ||
       prevMultiple( 9UL, 4UL ) != 8UL ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: prevMultiple() with 'unsigned int' failed\n"
          << " Details:\n"
          << "   prevMultiple( 7, 4 ) = " << prevMultiple( 7UL, 4UL ) << " (expected 4)\n"
          << "   prevMultiple( 8, 4 ) = " << prevMultiple( 8UL, 4UL ) << " (expected 8)\n"
          << "   prevMultiple( 9, 4 ) = " << prevMultiple( 9UL, 4UL ) << " (expected 8)\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************

} // namespace shims

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
   std::cout << "   Running mathematical shims operation test..." << std::endl;

   try
   {
      RUN_SHIMS_OPERATION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during mathematical shims operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
