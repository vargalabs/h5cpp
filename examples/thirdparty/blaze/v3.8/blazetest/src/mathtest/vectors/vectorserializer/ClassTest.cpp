//=================================================================================================
/*!
//  \file src/mathtest/vectors/vectorserializer/ClassTest.cpp
//  \brief Source file for the VectorSerializer class test
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
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/serialization/MatrixSerializer.h>
#include <blaze/math/serialization/VectorSerializer.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/vectors/vectorserializer/ClassTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace vectors {

namespace vectorserializer {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the VectorSerializer class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testEmptyVectors();
   testRandomVectors();
   testFailures();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Serialization test with empty vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs serialization tests with empty vectors. In case an error is
// detected, i.e. in case the destination vector is not empty after deserialization, a
// \a std::runtime_error exception is thrown.
*/
void ClassTest::testEmptyVectors()
{
   test_ = "Empty vectors";

   {
      blaze::DynamicVector<int> src;

      runDynamicVectorTests   ( src );
      runCompressedVectorTests( src );
   }

   {
      blaze::DynamicVector<int> src;

      runDynamicVectorTests   ( src );
      runCompressedVectorTests( src );
   }

   {
      blaze::CompressedVector<int> src;

      runDynamicVectorTests   ( src );
      runCompressedVectorTests( src );
   }

   {
      blaze::CompressedVector<int> src;

      runDynamicVectorTests   ( src );
      runCompressedVectorTests( src );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Serialization test with randomly initialized vectors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs serialization tests with randomly initialized vectors. In
// case an error is detected, i.e. in case a vector cannot be reconstituted from file,
// a \a std::runtime_error exception is thrown.
*/
void ClassTest::testRandomVectors()
{
   test_ = "Randomly initialized vectors";


   //=====================================================================================
   // StaticVector source
   //=====================================================================================

   {
      blaze::StaticVector<int,13UL> src;
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::StaticVector<unsigned int,13UL> src;
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::StaticVector<blaze::complex<float>,13UL> src;
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::StaticVector<blaze::StaticVector<double,3UL>,13UL> src;
      randomize( src );
      runAllTests<13UL>( src );
   }


   //=====================================================================================
   // DynamicVector source
   //=====================================================================================

   {
      blaze::DynamicVector<int> src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::DynamicVector<unsigned int> src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::DynamicVector< blaze::complex<float> > src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::DynamicVector< blaze::StaticVector<double,3UL> > src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }


   //=====================================================================================
   // CompressedVector source
   //=====================================================================================

   {
      blaze::CompressedVector<int> src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::CompressedVector<unsigned int> src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::CompressedVector< blaze::complex<float> > src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }

   {
      blaze::CompressedVector< blaze::StaticVector<double,3UL> > src( 13UL );
      randomize( src );
      runAllTests<13UL>( src );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of failing serialization attempts.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs tests with failing serialization attempts. In case no error is
// detected, i.e. in case the test is failing, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testFailures()
{
   test_ = "Serialization failures";

   try {
      blaze::DynamicMatrix<int> src( 5UL, 1UL );
      blaze::DynamicVector<int> dst;

      randomize( src );

      blaze::Archive<std::stringstream> archive;
      archive << src;
      archive >> dst;

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Type difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}

   try {
      blaze::DynamicVector<int> src( 5UL );
      blaze::StaticVector<int,3UL> dst;

      randomize( src );
      runTest( src, dst );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Size difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}

   try {
      blaze::DynamicVector<int> src( 5UL );
      blaze::DynamicVector<float> dst;

      randomize( src );
      runTest( src, dst );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Element type difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}

   try {
      blaze::DynamicVector<short> src( 5UL );
      blaze::DynamicVector<long int> dst;

      randomize( src );
      runTest( src, dst );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Element size difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}
}
//*************************************************************************************************

} // namespace vectorserializer

} // namespace vectors

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
   std::cout << "   Running VectorSerializer class test..." << std::endl;

   try
   {
      RUN_VECTORSERIALIZER_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during VectorSerializer class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
