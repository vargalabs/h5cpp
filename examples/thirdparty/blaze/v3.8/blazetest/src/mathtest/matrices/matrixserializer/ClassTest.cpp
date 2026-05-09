//=================================================================================================
/*!
//  \file src/mathtest/matrices/matrixserializer/ClassTest.cpp
//  \brief Source file for the MatrixSerializer class test
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
#include <blaze/math/DynamicVector.h>
#include <blaze/math/serialization/MatrixSerializer.h>
#include <blaze/math/serialization/VectorSerializer.h>
#include <blaze/math/StaticVector.h>
#include <blaze/util/Complex.h>
#include <blazetest/mathtest/matrices/matrixserializer/ClassTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace matrixserializer {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the MatrixSerializer class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testEmptyMatrices();
   testRandomMatrices();
   testFailures();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Serialization test with empty matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs serialization tests with empty matrices. In case an error is
// detected, i.e. in case the destination matrix is not empty after deserialization, a
// \a std::runtime_error exception is thrown.
*/
void ClassTest::testEmptyMatrices()
{
   test_ = "Empty matrices";

   {
      blaze::DynamicMatrix<int> src;

      runDynamicMatrixTests   ( src );
      runCompressedMatrixTests( src );
   }

   {
      blaze::DynamicMatrix<int> src;

      runDynamicMatrixTests   ( src );
      runCompressedMatrixTests( src );
   }

   {
      blaze::CompressedMatrix<int> src;

      runDynamicMatrixTests   ( src );
      runCompressedMatrixTests( src );
   }

   {
      blaze::CompressedMatrix<int> src;

      runDynamicMatrixTests   ( src );
      runCompressedMatrixTests( src );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Serialization test with randomly initialized matrices.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs serialization tests with randomly initialized matrices. In
// case an error is detected, i.e. in case a matrix cannot be reconstituted from file,
// a \a std::runtime_error exception is thrown.
*/
void ClassTest::testRandomMatrices()
{
   test_ = "Randomly initialized matrices";


   //=====================================================================================
   // Row-major StaticMatrix source
   //=====================================================================================

   {
      blaze::StaticMatrix<int,7UL,13UL,blaze::rowMajor> src;
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::StaticMatrix<unsigned int,7UL,13UL,blaze::rowMajor> src;
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::StaticMatrix<blaze::complex<float>,13UL,7UL,blaze::rowMajor> src;
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::StaticMatrix<blaze::StaticVector<double,3UL>,13UL,7UL,blaze::rowMajor> src;
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }


   //=====================================================================================
   // Column-major StaticMatrix source
   //=====================================================================================

   {
      blaze::StaticMatrix<int,7UL,13UL,blaze::columnMajor> src;
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::StaticMatrix<unsigned int,7UL,13UL,blaze::columnMajor> src;
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::StaticMatrix<blaze::complex<float>,13UL,7UL,blaze::columnMajor> src;
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::StaticMatrix<blaze::StaticVector<double,3UL>,13UL,7UL,blaze::columnMajor> src;
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }


   //=====================================================================================
   // Row-major DynamicMatrix source
   //=====================================================================================

   {
      blaze::DynamicMatrix<int,blaze::rowMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::DynamicMatrix<unsigned int,blaze::rowMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::DynamicMatrix<blaze::complex<float>,blaze::rowMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::DynamicMatrix<blaze::StaticVector<double,3UL>,blaze::rowMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }


   //=====================================================================================
   // Column-major DynamicMatrix source
   //=====================================================================================

   {
      blaze::DynamicMatrix<int,blaze::columnMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::DynamicMatrix<unsigned int,blaze::columnMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::DynamicMatrix<blaze::complex<float>,blaze::columnMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::DynamicMatrix<blaze::StaticVector<double,3UL>,blaze::columnMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }


   //=====================================================================================
   // Row-major CompressedMatrix source
   //=====================================================================================

   {
      blaze::CompressedMatrix<int,blaze::rowMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::CompressedMatrix<unsigned int,blaze::rowMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::CompressedMatrix<blaze::complex<float>,blaze::rowMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::CompressedMatrix<blaze::StaticVector<double,3UL>,blaze::rowMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }


   //=====================================================================================
   // Column-major CompressedMatrix source
   //=====================================================================================

   {
      blaze::CompressedMatrix<int,blaze::columnMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::CompressedMatrix<unsigned int,blaze::columnMajor> src( 7UL, 13UL );
      randomize( src );
      runAllTests<7UL,13UL>( src );
   }

   {
      blaze::CompressedMatrix<blaze::complex<float>,blaze::columnMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
   }

   {
      blaze::CompressedMatrix<blaze::StaticVector<double,3UL>,blaze::columnMajor> src( 13UL, 7UL );
      randomize( src );
      runAllTests<13UL,7UL>( src );
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
      blaze::DynamicVector<int> src( 10UL );
      blaze::DynamicMatrix<int> dst;

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
      blaze::DynamicMatrix<int> src( 4UL, 4UL );
      blaze::StaticMatrix<int,3UL,4UL> dst;

      randomize( src );
      runTest( src, dst );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Row difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}

   try {
      blaze::DynamicMatrix<int> src( 3UL, 5UL );
      blaze::StaticMatrix<int,3UL,4UL> dst;

      randomize( src );
      runTest( src, dst );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Column difference succeeded\n"
          << " Details:\n"
          << "   Source:\n" << src << "\n"
          << "   Destination:\n" << dst << "\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::runtime_error& )
   {}

   try {
      blaze::DynamicMatrix<int> src( 5UL, 4UL );
      blaze::DynamicMatrix<float> dst;

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
      blaze::DynamicMatrix<short> src( 5UL, 4UL );
      blaze::DynamicMatrix<long int> dst;

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

} // namespace matrixserializer

} // namespace matrices

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
   std::cout << "   Running MatrixSerializer class test..." << std::endl;

   try
   {
      RUN_MATRIXSERIALIZER_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during MatrixSerializer class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
