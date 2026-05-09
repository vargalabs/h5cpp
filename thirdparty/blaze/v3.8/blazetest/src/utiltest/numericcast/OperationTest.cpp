//=================================================================================================
/*!
//  \file src/utiltest/numericcast/OperationTest.cpp
//  \brief Source file for the numeric cast operation test
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

#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <blaze/util/NumericCast.h>
#include <blazetest/utiltest/numericcast/OperationTest.h>


namespace blazetest {

namespace utiltest {

namespace numericcast {

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
   // Integral/integral conversions
   testIntToInt();
   testIntToUInt();
   testUIntToInt();
   testIntToLong();
   testLongToInt();
   testULongToUInt();

   // Integral/floating point conversions
   testIntToFloat();
   testFloatToInt();
   testFloatToUInt();
   testIntToDouble();
   testDoubleToInt();

   // Floating point/floating point conversions
   testFloatToDouble();
   testDoubleToFloat();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of a numeric cast from \c int to \c int.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c int to \c int. In case an error is detected,
// a compilation error is created.
*/
void OperationTest::testIntToInt()
{
   {
      test_ = "Successful conversion from 'int' to 'int'";

      const int a( 2 );
      const int b( blaze::numeric_cast<int>( a ) );

      if( b != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = 2\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c int to \c unsigned int.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c int to \c unsigned int. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testIntToUInt()
{
   {
      test_ = "Successful conversion from 'int' to 'unsigned int'";

      const int a( 2 );
      const unsigned int b( blaze::numeric_cast<unsigned int>( a ) );

      if( b != 2U ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = 2\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Failing conversion from 'int' to 'unsigned int' (underflow)";

      const int a( -1 );

      try {
         const unsigned int b( blaze::numeric_cast<unsigned int>( a ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Successful numeric cast detected\n"
             << " Details:\n"
             << "   Source = " << a << "\n"
             << "   Result = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::underflow_error& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c unsigned int to \c int.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c unsigned int to \c int. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testUIntToInt()
{
   {
      test_ = "Successful conversion from 'unsigned int' to 'int'";

      const unsigned int a( 2U );
      const int b( blaze::numeric_cast<int>( a ) );

      if( b != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = 2\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Failing conversion from 'unsigned int' to 'int' (overflow)";

      constexpr unsigned int max( std::numeric_limits<int>::max() );

      const unsigned int a( max + 1U );

      try {
         const int b( blaze::numeric_cast<int>( a ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Successful numeric cast detected\n"
             << " Details:\n"
             << "   Source = " << a << "\n"
             << "   Result = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::overflow_error& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c int to \c long.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c int to \c long. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIntToLong()
{
   {
      test_ = "Successful conversion from 'int' to 'long'";

      const int a( 2 );
      const long b( blaze::numeric_cast<long>( a ) );

      if( b != 2L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = 2\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c long to \c int.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c long to \c int. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testLongToInt()
{
   {
      test_ = "Successful conversion from 'long' to 'int'";

      const long a( std::numeric_limits<int>::max() );
      const int b( blaze::numeric_cast<int>( a ) );

      if( b != std::numeric_limits<int>::max() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = " << std::numeric_limits<int>::max() << "\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Failing conversion from 'long' to 'int' (underflow)";

      constexpr long min( std::numeric_limits<int>::min() );

      const long a( min - 1L );

      try {
         const int b( blaze::numeric_cast<int>( a ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Successful numeric cast detected\n"
             << " Details:\n"
             << "   Source = " << a << "\n"
             << "   Result = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::underflow_error& ) {}
   }

   {
      test_ = "Failing conversion from 'long' to 'int' (overflow)";

      constexpr long max( std::numeric_limits<int>::max() );

      const long a( max + 1L );

      try {
         const int b( blaze::numeric_cast<int>( a ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Successful numeric cast detected\n"
             << " Details:\n"
             << "   Source = " << a << "\n"
             << "   Result = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::overflow_error& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c unsigned long to \c unsigned int.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c unsigned long to \c unsigned int. In case
// an error is detected, a compilation error is created.
*/
void OperationTest::testULongToUInt()
{
   {
      test_ = "Successful conversion from 'unsigned long' to 'unsigned int'";

      const unsigned long a( std::numeric_limits<unsigned int>::max() );
      const unsigned int b( blaze::numeric_cast<unsigned int>( a ) );

      if( b != std::numeric_limits<unsigned int>::max() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = " << std::numeric_limits<unsigned int>::max() << "\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Failing conversion from 'unsigned long' to 'unsigned int' (overflow)";

      constexpr unsigned long max( std::numeric_limits<unsigned int>::max() );

      const unsigned long a( max + 1UL );

      try {
         const unsigned int b( blaze::numeric_cast<unsigned int>( a ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Successful numeric cast detected\n"
             << " Details:\n"
             << "   Source = " << a << "\n"
             << "   Result = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::overflow_error& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c int to \c float.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c int to \c float. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIntToFloat()
{
   {
      test_ = "Successful conversion from 'int' to 'float'";

      const int a( 123 );
      const float b( blaze::numeric_cast<float>( a ) );

      if( b != 123.0F ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = 123\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c float to \c int.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c float to \c int. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testFloatToInt()
{
   {
      test_ = "Successful conversion from 'float' to 'int'";

      const float a( 123.4F );
      const int b( blaze::numeric_cast<int>( a ) );

      if( b != 123 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = 123\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Failing conversion from 'float' to 'int' (underflow)";

      const float a( -std::numeric_limits<float>::max() );

      try {
         const int b( blaze::numeric_cast<int>( a ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Successful numeric cast detected\n"
             << " Details:\n"
             << "   Source = " << a << "\n"
             << "   Result = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::underflow_error& ) {}
   }

   {
      test_ = "Failing conversion from 'float' to 'int' (overflow)";

      const float a( std::numeric_limits<float>::max() );

      try {
         const int b( blaze::numeric_cast<int>( a ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Successful numeric cast detected\n"
             << " Details:\n"
             << "   Source = " << a << "\n"
             << "   Result = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::overflow_error& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c float to \c unsigned int.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c float to \c unsigned int. In case an error
// is detected, a compilation error is created.
*/
void OperationTest::testFloatToUInt()
{
   {
      test_ = "Successful conversion from 'float' to 'unsigned int'";

      const float a( 123.4F );
      const unsigned int b( blaze::numeric_cast<unsigned int>( a ) );

      if( b != 123 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = 123\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Failing conversion from 'float' to 'unsigned int' (underflow)";

      const float a( -1.0F );

      try {
         const unsigned int b( blaze::numeric_cast<unsigned int>( a ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Successful numeric cast detected\n"
             << " Details:\n"
             << "   Source = " << a << "\n"
             << "   Result = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::underflow_error& ) {}
   }

   {
      test_ = "Failing conversion from 'float' to 'unsigned int' (overflow)";

      const float a( std::numeric_limits<float>::max() );

      try {
         const unsigned int b( blaze::numeric_cast<unsigned int>( a ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Successful numeric cast detected\n"
             << " Details:\n"
             << "   Source = " << a << "\n"
             << "   Result = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::overflow_error& ) {}
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c int to \c double.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c int to \c double. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testIntToDouble()
{
   {
      test_ = "Successful conversion from 'int' to 'double'";

      const int a( 123 );
      const double b( blaze::numeric_cast<double>( a ) );

      if( b != 123.0F ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = 123\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c double to \c int.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c double to \c int. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testDoubleToInt()
{
   {
      test_ = "Successful conversion from 'double' to 'int'";

      const double a( 123.4 );
      const int b( blaze::numeric_cast<int>( a ) );

      if( b != 123 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = 123\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c float to \c double.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c float to \c double. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testFloatToDouble()
{
   {
      test_ = "Successful conversion from 'float' to 'double'";

      const float a( std::numeric_limits<float>::max() );
      const double b( blaze::numeric_cast<double>( a ) );

      if( b != std::numeric_limits<float>::max() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = 123\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a numeric cast from \c double to \c float.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs several conversions from \c double to \c float. In case an error is
// detected, a compilation error is created.
*/
void OperationTest::testDoubleToFloat()
{
   {
      test_ = "Successful conversion from 'double' to 'float'";

      const double a( std::numeric_limits<float>::max() );
      const float b( blaze::numeric_cast<float>( a ) );

      if( b != std::numeric_limits<float>::max() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed numeric cast detected\n"
             << " Details:\n"
             << "   Expected result = " << std::numeric_limits<float>::max() << "\n"
             << "   Actual result   = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Failing conversion from 'double' to 'float' (underflow)";

      const double a( -std::numeric_limits<double>::max() );

      try {
         const float b( blaze::numeric_cast<float>( a ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Successful numeric cast detected\n"
             << " Details:\n"
             << "   Source = " << a << "\n"
             << "   Result = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::underflow_error& ) {}
   }

   {
      test_ = "Failing conversion from 'double' to 'float' (overflow)";

      const double a( std::numeric_limits<double>::max() );

      try {
         const float b( blaze::numeric_cast<float>( a ) );

         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Successful numeric cast detected\n"
             << " Details:\n"
             << "   Source = " << a << "\n"
             << "   Result = " << b << "\n";
         throw std::runtime_error( oss.str() );
      }
      catch( std::overflow_error& ) {}
   }
}
//*************************************************************************************************

} // namespace numericcast

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
   std::cout << "   Running numeric cast operation test..." << std::endl;

   try
   {
      RUN_NUMERICCAST_OPERATION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during numeric cast operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
