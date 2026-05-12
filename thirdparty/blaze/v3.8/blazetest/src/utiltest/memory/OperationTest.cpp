//=================================================================================================
/*!
//  \file src/utiltest/memory/OperationTest.cpp
//  \brief Source file for the memory operation test
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
#include <blaze/math/StaticVector.h>
#include <blaze/util/Memory.h>
#include <blaze/util/typetraits/AlignmentOf.h>
#include <blazetest/utiltest/memory/OperationTest.h>
#include <blazetest/utiltest/AlignedResource.h>
#include <blazetest/utiltest/ThrowingResource.h>


namespace blazetest {

namespace utiltest {

namespace memory {

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
   testBuiltinTypes();
   testClassTypes();
   testNullPointer();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST TYPE TRAITS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of aligned allocation and deallocation of built-in data types.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the aligned allocation and deallocation functionality in
// combination with built-in data types. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void OperationTest::testBuiltinTypes()
{
   // Built-in data type 'char'
   {
      test_ = "Built-in data types (char)";

      char* array = blaze::allocate<char>( number );

      const size_t alignment( blaze::AlignmentOf<char>::value );
      const size_t deviation( reinterpret_cast<size_t>( array ) % alignment );

      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation         : " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }

      blaze::deallocate( array );
   }

   // Built-in data type 'int'
   {
      test_ = "Built-in data types (int)";

      int* array = blaze::allocate<int>( number );

      const size_t alignment( blaze::AlignmentOf<int>::value );
      const size_t deviation( reinterpret_cast<size_t>( array ) % alignment );

      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation         : " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }

      blaze::deallocate( array );
   }

   // Built-in data type 'float'
   {
      test_ = "Built-in data types (float)";

      float* array = blaze::allocate<float>( number );

      const size_t alignment( blaze::AlignmentOf<float>::value );
      const size_t deviation( reinterpret_cast<size_t>( array ) % alignment );

      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation         : " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }

      blaze::deallocate( array );
   }

   // Built-in data type 'double'
   {
      test_ = "Built-in data types (double)";

      double* array = blaze::allocate<double>( number );

      const size_t alignment( blaze::AlignmentOf<double>::value );
      const size_t deviation( reinterpret_cast<size_t>( array ) % alignment );

      if( deviation != 0UL ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid alignment detected\n"
             << " Details:\n"
             << "   Expected alignment: " << alignment << "\n"
             << "   Deviation         : " << deviation << "\n";
         throw std::runtime_error( oss.str() );
      }

      blaze::deallocate( array );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of aligned allocation and deallocation of class types.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the aligned allocation and deallocation functionality in
// combination with class types. In case an error is detected, a \a std::runtime_error exception
// is thrown.
*/
void OperationTest::testClassTypes()
{
   // User-specific class type 'AlignedResource'
   {
      test_ = "User-specific class types (AlignedResource)";

      AlignedResource* array = blaze::allocate<AlignedResource>( number );

      const size_t alignment( blaze::AlignmentOf<AlignedResource>::value );

      for( size_t i=0UL; i<number; ++i )
      {
         const size_t deviation( reinterpret_cast<size_t>( &array[i] ) % alignment );

         if( deviation != 0UL ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid alignment at index " << i << " detected\n"
                << " Details:\n"
                << "   Expected alignment: " << alignment << "\n"
                << "   Deviation         : " << deviation << "\n";
            throw std::runtime_error( oss.str() );
         }

         if( array[i].getValue() != 7U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid value at index " << i << " detected\n"
                << " Details:\n"
                << "   Current value : " << array[i].getValue() << "\n"
                << "   Expected value: 7\n";
            throw std::runtime_error( oss.str() );
         }
      }

      blaze::deallocate( array );

      if( AlignedResource::getCount() != 0U ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of instances detected\n"
             << " Details:\n"
             << "   Current count : " << AlignedResource::getCount() << "\n"
             << "   Expected count: 0\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // User-specific class type 'ThrowingResource'
   {
      test_ = "User-specific class types (ThrowingResource)";

      try {
         blaze::allocate<ThrowingResource>( number );
      }
      catch( ... )
      {
         if( ThrowingResource::getCount() != 0U ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Invalid number of instances detected\n"
                << " Details:\n"
                << "   Current count : " << ThrowingResource::getCount() << "\n"
                << "   Expected count: 0\n";
            throw std::runtime_error( oss.str() );
         }
         return;
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of a deallocation of a \c nullptr.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the deallocation using a \c nullptr. In case an error is
// detected, a \a std::runtime_error exception is thrown.
*/
void OperationTest::testNullPointer()
{
   const int* const array( nullptr );

   blaze::deallocate( array );
}
//*************************************************************************************************

} // namespace memory

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
   std::cout << "   Running memory operation test..." << std::endl;

   try
   {
      RUN_MEMORY_OPERATION_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during memory operation test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
