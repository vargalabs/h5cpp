//=================================================================================================
/*!
//  \file src/mathtest/vectors/staticvector/ClassTest.cpp
//  \brief Source file for the StaticVector class test
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

#include <array>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/CustomVector.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/util/Complex.h>
#include <blaze/util/policies/Deallocate.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>
#include <blazetest/mathtest/vectors/staticvector/ClassTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace vectors {

namespace staticvector {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the StaticVector class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testAlignment< char           >( "char"           );
   testAlignment< signed char    >( "signed char"    );
   testAlignment< unsigned char  >( "unsigned char"  );
   testAlignment< wchar_t        >( "wchar_t"        );
   testAlignment< short          >( "short"          );
   testAlignment< unsigned short >( "unsigned short" );
   testAlignment< int            >( "int"            );
   testAlignment< unsigned int   >( "unsigned int"   );
   testAlignment< long           >( "long"           );
   testAlignment< unsigned long  >( "unsigned long"  );
   testAlignment< float          >( "float"          );
   testAlignment< double         >( "double"         );

   testAlignment< complex<char>           >( "complex<char>"           );
   testAlignment< complex<signed char>    >( "complex<signed char>"    );
   testAlignment< complex<unsigned char>  >( "complex<unsigned char>"  );
   testAlignment< complex<wchar_t>        >( "complex<wchar_t>"        );
   testAlignment< complex<short>          >( "complex<short>"          );
   testAlignment< complex<unsigned short> >( "complex<unsigned short>" );
   testAlignment< complex<int>            >( "complex<int>"            );
   testAlignment< complex<unsigned int>   >( "complex<unsigned int>"   );
   testAlignment< complex<float>          >( "complex<float>"          );
   testAlignment< complex<double>         >( "complex<double>"         );

   testConstructors();
   testAssignment();
   testAddAssign();
   testSubAssign();
   testMultAssign();
   testDivAssign();
   testCrossAssign();
   testScaling();
   testSubscript();
   testAt();
   testIterator();
   testNonZeros();
   testReset();
   testClear();
   testSwap();
   testIsDefault();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the StaticVector constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the StaticVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Default constructor
   //=====================================================================================

   {
      test_ = "StaticVector default constructor (size 0)";

      blaze::StaticVector<int,0UL,blaze::rowVector> vec;

      checkSize    ( vec, 0UL );
      checkCapacity( vec, 0UL );
      checkNonZeros( vec, 0UL );
   }

   {
      test_ = "StaticVector default constructor (size 5)";

      blaze::StaticVector<int,5UL,blaze::rowVector> vec;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 0UL );

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Homogeneous initialization
   //=====================================================================================

   {
      test_ = "StaticVector homogeneous initialization constructor";

      blaze::StaticVector<int,3UL,blaze::rowVector> vec( 2 );

      checkSize    ( vec, 3UL );
      checkCapacity( vec, 3UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // List initialization
   //=====================================================================================

   {
      test_ = "StaticVector initializer list constructor (complete list)";

      blaze::StaticVector<int,4UL,blaze::rowVector> vec{ 1, 2, 3, 4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector initializer list constructor (incomplete list)";

      blaze::StaticVector<int,4UL,blaze::rowVector> vec{ 1, 2 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 2UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Array initialization
   //=====================================================================================

   {
      test_ = "StaticVector dynamic array initialization constructor";

      std::unique_ptr<int[]> array( new int[2] );
      array[0] = 1;
      array[1] = 2;
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( 2UL, array.get() );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 2UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dynamic array initialization constructor";

      std::unique_ptr<int[]> array( new int[4] );
      array[0] = 1;
      array[1] = 2;
      array[2] = 3;
      array[3] = 4;
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( 4UL, array.get() );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector static array initialization constructor";

      int array[4] = { 1, 2, 3, 4 };
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( array );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector std::array initialization constructor";

      std::array<int,4UL> array{ 1, 2, 3, 4 };
      blaze::StaticVector<int,4UL,blaze::rowVector> vec( array );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy constructor
   //=====================================================================================

   {
      test_ = "StaticVector copy constructor";

      blaze::StaticVector<int,5UL,blaze::rowVector> vec1{ 1, 2, 3, 4, 5 };
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Dense vector constructor
   //=====================================================================================

   {
      test_ = "StaticVector dense vector constructor (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[2] = 3;
      vec1[3] = 4;
      vec1[4] = 5;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector constructor (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[2] = 3;
      vec1[3] = 4;
      vec1[4] = 5;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      test_ = "StaticVector dense vector constructor (non-fitting vector)";

      blaze::DynamicVector<int,blaze::rowVector> vec1{ 1, 2, 3, 4, 5 };
      blaze::StaticVector<int,4UL,blaze::rowVector> vec2( vec1 );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Construction succeeded\n"
          << " Details:\n"
          << "   Result:\n" << vec2 << "\n"
          << "   Given vector:\n( 1 2 3 4 5 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& )
   {}


   //=====================================================================================
   // Sparse vector constructor
   //=====================================================================================

   {
      test_ = "StaticVector sparse vector constructor";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 5UL, 3UL );
      vec1[0] = 1;
      vec1[2] = 3;
      vec1[4] = 5;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2( vec1 );

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 0 || vec2[2] != 3 || vec2[3] != 0 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 0 3 0 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   try {
      test_ = "StaticVector sparse vector constructor (non-fitting vector)";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 5UL, 3UL );
      vec1[0] = 1;
      vec1[2] = 3;
      vec1[4] = 5;
      blaze::StaticVector<int,4UL,blaze::rowVector> vec2( vec1 );

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Construction succeeded\n"
          << " Details:\n"
          << "   Result:\n" << vec2 << "\n"
          << "   Given vector:\n( 1 0 3 0 5 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::invalid_argument& )
   {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the StaticVector assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the StaticVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Homogeneous assignment
   //=====================================================================================

   {
      test_ = "StaticVector homogeneous assignment";

      blaze::StaticVector<int,5UL,blaze::rowVector> vec;
      vec = 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 5UL );

      if( vec[0] != 2 || vec[1] != 2 || vec[2] != 2 || vec[3] != 2 || vec[4] != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // List assignment
   //=====================================================================================

   {
      test_ = "StaticVector initializer list assignment (complete list)";

      blaze::StaticVector<int,4UL,blaze::rowVector> vec;
      vec = { 1, 2, 3, 4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector initializer list assignment (incomplete list)";

      blaze::StaticVector<int,4UL,blaze::rowVector> vec;
      vec = { 1, 2 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 2UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Array assignment
   //=====================================================================================

   {
      test_ = "StaticVector static array assignment";

      int array[4] = { 1, 2 };
      blaze::StaticVector<int,4UL,blaze::rowVector> vec;
      vec = array;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 2UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector static array assignment";

      int array[4] = { 1, 2, 3, 4 };
      blaze::StaticVector<int,4UL,blaze::rowVector> vec;
      vec = array;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector std::array assignment";

      std::array<int,4UL> array{ 1, 2, 3, 4 };
      blaze::StaticVector<int,4UL,blaze::rowVector> vec;
      vec = array;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Copy assignment
   //=====================================================================================

   {
      test_ = "StaticVector copy assignment";

      blaze::StaticVector<int,5UL,blaze::rowVector> vec1{ 1, 2, 3, 4, 5 };
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector copy assignment stress test";

      using RandomVectorType = blaze::StaticVector<int,10UL,blaze::rowVector>;

      blaze::StaticVector<int,10UL,blaze::rowVector> vec1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const RandomVectorType vec2( blaze::rand<RandomVectorType>( min, max ) );

         vec1 = vec2;

         if( vec1 != vec2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << vec1 << "\n"
                << "   Expected result:\n" << vec2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Dense vector assignment
   //=====================================================================================

   {
      test_ = "StaticVector dense vector assignment (mixed type)";

      blaze::StaticVector<short,5UL,blaze::rowVector> vec1{ 1, 2, 3, 4, 5 };
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[2] = 3;
      vec1[3] = 4;
      vec1[4] = 5;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1[0] = 1;
      vec1[1] = 2;
      vec1[2] = 3;
      vec1[3] = 4;
      vec1[4] = 5;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 5UL );

      if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 || vec2[4] != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 2 3 4 5 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector copy assignment stress test";

      using RandomVectorType = blaze::DynamicVector<int,blaze::rowVector>;

      blaze::StaticVector<int,10UL,blaze::rowVector> vec1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const RandomVectorType vec2( blaze::rand<RandomVectorType>( 10UL, min, max ) );

         vec1 = vec2;

         if( vec1 != vec2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << vec1 << "\n"
                << "   Expected result:\n" << vec2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Sparse vector assignment
   //=====================================================================================

   {
      test_ = "StaticVector sparse vector assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 5UL );
      vec1[0] = 1;
      vec1[2] = 2;
      vec1[3] = 3;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2;
      vec2 = vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 1 || vec2[1] != 0 || vec2[2] != 2 || vec2[3] != 3 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 0 2 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector copy assignment stress test";

      using RandomVectorType = blaze::CompressedVector<int,blaze::rowVector>;

      blaze::StaticVector<int,10UL,blaze::rowVector> vec1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const RandomVectorType vec2( blaze::rand<RandomVectorType>( 10UL, min, max ) );

         vec1 = vec2;

         if( vec1 != vec2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << vec1 << "\n"
                << "   Expected result:\n" << vec2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the StaticVector addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the StaticVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // Dense vector addition assignment
   //=====================================================================================

   {
      test_ = "StaticVector dense vector addition assignment (mixed type)";

      blaze::StaticVector<short,5UL,blaze::rowVector> vec1{ 1, 0, -2, 3, 0 };
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4, 2, -6, 7 };

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1[0] =  1;
      vec1[1] =  0;
      vec1[2] = -2;
      vec1[3] =  3;
      vec1[4] =  0;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4, 2, -6, 7 };

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1[0] =  1;
      vec1[1] =  0;
      vec1[2] = -2;
      vec1[3] =  3;
      vec1[4] =  0;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4,  2, -6, 7 };

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "DynamicVector sparse vector addition assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 5UL, 3UL );
      vec1[0] =  1;
      vec1[2] = -2;
      vec1[3] =  3;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4, 2, -6, 7 };

      vec2 += vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the StaticVector subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the StaticVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // Dense vector subtraction assignment
   //=====================================================================================

   {
      test_ = "StaticVector dense vector subtraction assignment (mixed type)";

      blaze::StaticVector<short,5UL,blaze::rowVector> vec1{ -1, 0, 2, -3, 0 };
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4, 2, -6, 7 };

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1[0] = -1;
      vec1[1] =  0;
      vec1[2] =  2;
      vec1[3] = -3;
      vec1[4] =  0;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4, 2, -6, 7 };

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1[0] = -1;
      vec1[1] =  0;
      vec1[2] =  2;
      vec1[3] = -3;
      vec1[4] =  0;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4, 2, -6, 7 };

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector addition assignment
   //=====================================================================================

   {
      test_ = "StaticVector sparse vector subtraction assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 5UL, 3UL );
      vec1[0] = -1;
      vec1[2] =  2;
      vec1[3] = -3;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4, 2, -6, 7 };

      vec2 -= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 4UL );

      if( vec2[0] != 1 || vec2[1] != 4 || vec2[2] != 0 || vec2[3] != -3 || vec2[4] != 7 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 1 4 0 -3 7 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the StaticVector multiplication assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the multiplication assignment operators of the StaticVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testMultAssign()
{
   //=====================================================================================
   // Dense vector multiplication assignment
   //=====================================================================================

   {
      test_ = "StaticVector dense vector multiplication assignment (mixed type)";

      blaze::StaticVector<short,5UL,blaze::rowVector> vec1{ 1, 0, -2, 3, 0 };
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4,  2, -6, 7 };

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 2UL );

      if( vec2[0] != 0 || vec2[1] != 0 || vec2[2] != -4 || vec2[3] != -18 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 0 -4 -18 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector multiplication assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1[0] =  1;
      vec1[1] =  0;
      vec1[2] = -2;
      vec1[3] =  3;
      vec1[4] =  0;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4,  2, -6, 7 };

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 2UL );

      if( vec2[0] != 0 || vec2[1] != 0 || vec2[2] != -4 || vec2[3] != -18 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 0 -4 -18 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector multiplication assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1[0] =  1;
      vec1[1] =  0;
      vec1[2] = -2;
      vec1[3] =  3;
      vec1[4] =  0;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4,  2, -6, 7 };

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 2UL );

      if( vec2[0] != 0 || vec2[1] != 0 || vec2[2] != -4 || vec2[3] != -18 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 0 -4 -18 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector multiplication assignment
   //=====================================================================================

   {
      test_ = "StaticVector sparse vector multiplication assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 5UL, 3UL );
      vec1[0] =  1;
      vec1[2] = -2;
      vec1[3] =  3;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 0, 4, 2, -6, 7 };

      vec2 *= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 2UL );

      if( vec2[0] != 0 || vec2[1] != 0 || vec2[2] != -4 || vec2[3] != -18 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 0 -4 -18 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the StaticVector division assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the division assignment operators of the StaticVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testDivAssign()
{
   //=====================================================================================
   // Dense vector division assignment
   //=====================================================================================

   {
      test_ = "StaticVector dense vector division assignment (mixed type)";

      blaze::StaticVector<short,5UL,blaze::rowVector> vec1{ 1, 2, -3, 4, 1 };
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 2, 0, -3, 8, 0 };

      vec2 /= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 2 || vec2[1] != 0 || vec2[2] != 1 || vec2[3] != 2 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 0 1 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector division assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 5UL, 16UL );
      vec1[0] =  1;
      vec1[1] =  2;
      vec1[2] = -3;
      vec1[3] =  4;
      vec1[4] =  1;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 2, 0, -3, 8, 0 };

      vec2 /= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 2 || vec2[1] != 0 || vec2[2] != 1 || vec2[3] != 2 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 0 1 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector division assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[6] );
      UnalignedUnpadded vec1( memory.get()+1UL, 5UL );
      vec1[0] =  1;
      vec1[1] =  2;
      vec1[2] = -3;
      vec1[3] =  4;
      vec1[4] =  1;
      blaze::StaticVector<int,5UL,blaze::rowVector> vec2{ 2, 0, -3, 8, 0 };

      vec2 /= vec1;

      checkSize    ( vec2, 5UL );
      checkCapacity( vec2, 5UL );
      checkNonZeros( vec2, 3UL );

      if( vec2[0] != 2 || vec2[1] != 0 || vec2[2] != 1 || vec2[3] != 2 || vec2[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 2 0 1 2 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the StaticVector cross product assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the cross product assignment operators of the StaticVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testCrossAssign()
{
   //=====================================================================================
   // Dense vector cross product assignment
   //=====================================================================================

   {
      test_ = "StaticVector dense vector cross product assignment (mixed type)";

      blaze::StaticVector<short,3UL,blaze::rowVector> vec1{ 1, 0, -2 };
      blaze::StaticVector<int,3UL,blaze::rowVector> vec2{ 2, 0, -1 };

      vec2 %= vec1;

      checkSize    ( vec2, 3UL );
      checkCapacity( vec2, 3UL );
      checkNonZeros( vec2, 1UL );

      if( vec2[0] != 0 || vec2[1] != 3 || vec2[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector cross product assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowVector;

      using AlignedPadded = blaze::CustomVector<int,aligned,padded,rowVector>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 16UL ) );
      AlignedPadded vec1( memory.get(), 3UL, 16UL );
      vec1[0] =  1;
      vec1[1] =  0;
      vec1[2] = -2;
      blaze::StaticVector<int,3UL,blaze::rowVector> vec2{ 2, 0, -1 };

      vec2 %= vec1;

      checkSize    ( vec2, 3UL );
      checkCapacity( vec2, 3UL );
      checkNonZeros( vec2, 1UL );

      if( vec2[0] != 0 || vec2[1] != 3 || vec2[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector dense vector cross product assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowVector;

      using UnalignedUnpadded = blaze::CustomVector<int,unaligned,unpadded,rowVector>;
      std::unique_ptr<int[]> memory( new int[4] );
      UnalignedUnpadded vec1( memory.get()+1UL, 3UL );
      vec1[0] =  1;
      vec1[1] =  0;
      vec1[2] = -2;
      blaze::StaticVector<int,3UL,blaze::rowVector> vec2{ 2, 0, -1 };

      vec2 %= vec1;

      checkSize    ( vec2, 3UL );
      checkCapacity( vec2, 3UL );
      checkNonZeros( vec2, 1UL );

      if( vec2[0] != 0 || vec2[1] != 3 || vec2[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Sparse vector cross product assignment
   //=====================================================================================

   {
      test_ = "StaticVector sparse vector cross product assignment";

      blaze::CompressedVector<int,blaze::rowVector> vec1( 3UL, 2UL );
      vec1[0] =  1;
      vec1[2] = -2;
      blaze::StaticVector<int,3UL,blaze::rowVector> vec2{ 2, 0, -1 };

      vec2 %= vec1;

      checkSize    ( vec2, 3UL );
      checkCapacity( vec2, 3UL );
      checkNonZeros( vec2, 1UL );

      if( vec2[0] != 0 || vec2[1] != 3 || vec2[2] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Cross product assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << vec2 << "\n"
             << "   Expected result:\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of all StaticVector (self-)scaling operations.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all available ways to scale an instance of the StaticVector
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testScaling()
{
   //=====================================================================================
   // Self-scaling (v*=s)
   //=====================================================================================

   {
      test_ = "StaticVector self-scaling (v*=s)";

      blaze::StaticVector<int,5UL,blaze::rowVector> vec{ 1, 0, -2, 3, 0 };

      vec *= 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 0 || vec[2] != -4 || vec[3] != 6 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 0 -4 6 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v*s)
   //=====================================================================================

   {
      test_ = "StaticVector self-scaling (v=v*s)";

      blaze::StaticVector<int,5UL,blaze::rowVector> vec{ 1, 0, -2, 3, 0 };

      vec = vec * 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 0 || vec[2] != -4 || vec[3] != 6 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 0 -4 6 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=s*v)
   //=====================================================================================

   {
      test_ = "StaticVector self-scaling (v=s*v)";

      blaze::StaticVector<int,5UL,blaze::rowVector> vec{ 1, 0, -2, 3, 0 };

      vec = 2 * vec;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 2 || vec[1] != 0 || vec[2] != -4 || vec[3] != 6 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 0 -4 6 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v/=s)
   //=====================================================================================

   {
      test_ = "StaticVector self-scaling (v/=s)";

      blaze::StaticVector<int,5UL,blaze::rowVector> vec{ 2, 0, -4, 6, 0 };

      vec /= 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 1 || vec[1] != 0 || vec[2] != -2 || vec[3] != 3 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 -2 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Self-scaling (v=v/s)
   //=====================================================================================

   {
      test_ = "StaticVector self-scaling (v=v/s)";

      blaze::StaticVector<int,5UL,blaze::rowVector> vec{ 2, 0, -4, 6, 0 };

      vec = vec / 2;

      checkSize    ( vec, 5UL );
      checkCapacity( vec, 5UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 1 || vec[1] != 0 || vec[2] != -2 || vec[3] != 3 || vec[4] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed self-scaling operation\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 0 -2 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // StaticVector::scale()
   //=====================================================================================

   {
      test_ = "StaticVector::scale() (int)";

      // Initialization check
      blaze::StaticVector<int,4UL,blaze::rowVector> vec{ 1, 2, 3, 4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Integral scaling of the vector
      vec.scale( 2 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 2 || vec[1] != 4 || vec[2] != 6 || vec[3] != 8 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 2 4 6 8 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Floating point scaling of the vector
      vec.scale( 0.5 );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "StaticVector::scale() (complex)";

      using blaze::complex;

      blaze::StaticVector<complex<float>,2UL,blaze::rowVector> vec;
      vec[0] = complex<float>( 1.0F, 0.0F );
      vec[1] = complex<float>( 2.0F, 0.0F );
      vec.scale( complex<float>( 3.0F, 0.0F ) );

      checkSize    ( vec, 2UL );
      checkCapacity( vec, 2UL );
      checkNonZeros( vec, 2UL );

      if( vec[0] != complex<float>( 3.0F, 0.0F ) || vec[1] != complex<float>( 6.0F, 0.0F ) ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scale operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( (3,0) (6,0) )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the StaticVector subscript operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the subscript operator
// of the StaticVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testSubscript()
{
   test_ = "StaticVector::operator[]";

   // Assignment to the element at index 2
   blaze::StaticVector<int,5UL,blaze::rowVector> vec;
   vec[2] = 1;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 1UL );

   if( vec[2] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 4
   vec[4] = 2;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 2UL );

   if( vec[2] != 1 || vec[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 3
   vec[3] = 3;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 3UL );

   if( vec[2] != 1 || vec[3] != 3 || vec[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 3 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 0
   vec[0] = 4;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 4UL );

   if( vec[0] != 4 || vec[2] != 1 || vec[3] != 3 || vec[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 1 3 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Addition assignment to the element at index 2
   vec[2] += vec[3];

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 4UL );

   if( vec[0] != 4 || vec[2] != 4 || vec[3] != 3 || vec[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 4 3 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Subtraction assignment to the element at index 1
   vec[1] -= vec[4];

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != 3 || vec[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 4 3 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Multiplication assignment to the element at index 3
   vec[3] *= -3;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec[0] != 4 || vec[1] != -2 || vec[2] != 4 || vec[3] != -9 || vec[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 4 -9 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Division assignment to the element at index 2
   vec[2] /= 2;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec[0] != 4 || vec[1] != -2 || vec[2] != 2 || vec[3] != -9 || vec[4] != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Subscript operator failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 2 -9 2 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the StaticVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the StaticVector class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testAt()
{
   test_ = "StaticVector::at()";

   // Assignment to the element at index 2
   blaze::StaticVector<int,5UL,blaze::rowVector> vec;
   vec.at(2) = 1;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 1UL );

   if( vec.at(2) != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 4
   vec.at(4) = 2;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 2UL );

   if( vec.at(2) != 1 || vec.at(4) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 0 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 3
   vec.at(3) = 3;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 3UL );

   if( vec.at(2) != 1 || vec.at(3) != 3 || vec.at(4) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 0 0 1 3 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Assignment to the element at index 0
   vec.at(0) = 4;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 4UL );

   if( vec.at(0) != 4 || vec.at(2) != 1 || vec.at(3) != 3 || vec.at(4) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 1 3 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Addition assignment to the element at index 2
   vec.at(2) += vec.at(3);

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 4UL );

   if( vec.at(0) != 4 || vec.at(2) != 4 || vec.at(3) != 3 || vec.at(4) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 0 4 3 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Subtraction assignment to the element at index 1
   vec.at(1) -= vec.at(4);

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec.at(0) != 4 || vec.at(1) != -2 || vec.at(2) != 4 || vec.at(3) != 3 || vec.at(4) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 4 3 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Multiplication assignment to the element at index 3
   vec.at(3) *= -3;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec.at(0) != 4 || vec.at(1) != -2 || vec.at(2) != 4 || vec.at(3) != -9 || vec.at(4) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 4 -9 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Division assignment to the element at index 2
   vec.at(2) /= 2;

   checkSize    ( vec, 5UL );
   checkCapacity( vec, 5UL );
   checkNonZeros( vec, 5UL );

   if( vec.at(0) != 4 || vec.at(1) != -2 || vec.at(2) != 2 || vec.at(3) != -9 || vec.at(4) != 2 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Access via at() function failed\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 2 -9 2 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Attempt to assign to the element at index 5
   try {
      vec.at(5) = 2;

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Result:\n" << vec << "\n"
          << "   Expected result:\n( 4 -2 2 -9 2 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the StaticVector iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the StaticVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   using VectorType    = blaze::StaticVector<int,4UL>;
   using Iterator      = VectorType::Iterator;
   using ConstIterator = VectorType::ConstIterator;

   VectorType vec{ 1, 0, -2, -3 };

   // Testing the Iterator default constructor
   {
      test_ = "Iterator default constructor";

      Iterator it{};

      if( it != Iterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing the ConstIterator default constructor
   {
      test_ = "ConstIterator default constructor";

      ConstIterator it{};

      if( it != ConstIterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing conversion from Iterator to ConstIterator
   {
      test_ = "Iterator/ConstIterator conversion";

      ConstIterator it( begin( vec ) );

      if( it == end( vec ) || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator conversion detected\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via Iterator (end-begin)
   {
      test_ = "Iterator subtraction (end-begin)";

      const ptrdiff_t number( end( vec ) - begin( vec ) );

      if( number != 4L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via Iterator (begin-end)
   {
      test_ = "Iterator subtraction (begin-end)";

      const ptrdiff_t number( begin( vec ) - end( vec ) );

      if( number != -4L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: -4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via ConstIterator (end-begin)
   {
      test_ = "ConstIterator subtraction (end-begin)";

      const ptrdiff_t number( cend( vec ) - cbegin( vec ) );

      if( number != 4L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements via ConstIterator (begin-end)
   {
      test_ = "ConstIterator subtraction (begin-end)";

      const ptrdiff_t number( cbegin( vec ) - cend( vec ) );

      if( number != -4L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: -4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      ConstIterator it ( cbegin( vec ) );
      ConstIterator end( cend( vec ) );

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      --it;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it++;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it--;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it += 2UL;

      if( it == end || *it != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator addition assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it -= 2UL;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator subtraction assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it + 3UL;

      if( it == end || *it != -3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar addition failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it - 3UL;

      if( it == end || *it != 1 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar subtraction failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = 4UL + it;

      if( it != end ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scalar/iterator addition failed\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing assignment via Iterator
   {
      test_ = "Assignment via Iterator";

      int value = 6;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it = value++;
      }

      if( vec[0] != 6 || vec[1] != 7 || vec[2] != 8 || vec[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing addition assignment via Iterator
   {
      test_ = "Addition assignment via Iterator";

      int value = 2;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it += value++;
      }

      if( vec[0] != 8 || vec[1] != 10 || vec[2] != 12 || vec[3] != 14 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 8 10 12 14 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing subtraction assignment via Iterator
   {
      test_ = "Subtraction assignment via Iterator";

      int value = 2;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it -= value++;
      }

      if( vec[0] != 6 || vec[1] != 7 || vec[2] != 8 || vec[3] != 9 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 6 7 8 9 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing multiplication assignment via Iterator
   {
      test_ = "Multiplication assignment via Iterator";

      int value = 1;

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it *= value++;
      }

      if( vec[0] != 6 || vec[1] != 14 || vec[2] != 24 || vec[3] != 36 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Multiplication assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 6 14 24 36 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing division assignment via Iterator
   {
      test_ = "Division assignment via Iterator";

      for( Iterator it=begin( vec ); it!=end( vec ); ++it ) {
         *it /= 2;
      }

      if( vec[0] != 3 || vec[1] != 7 || vec[2] != 12 || vec[3] != 18 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Division assignment via iterator failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 3 7 12 18 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the StaticVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the StaticVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   test_ = "StaticVector::nonZeros()";

   {
      blaze::StaticVector<int,4UL,blaze::rowVector> vec;

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 0UL );

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      blaze::StaticVector<int,4UL,blaze::rowVector> vec{ 1, 2, 0, 3 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 3 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 0 3 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c reset() member function of the StaticVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c reset() member function of the StaticVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testReset()
{
   using blaze::reset;

   test_ = "StaticVector::reset()";

   // Resetting a default constructed vector
   {
      blaze::StaticVector<int,4UL,blaze::rowVector> vec;

      reset( vec );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 0UL );

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Resetting an initialized vector
   {
      // Initialization check
      blaze::StaticVector<int,4UL,blaze::rowVector> vec{ 1, 2, 3, 4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting a single element
      reset( vec[2] );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Resetting the vector
      reset( vec );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 0UL );

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Reset operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c clear() function with the StaticVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c clear() function with the StaticVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testClear()
{
   using blaze::clear;

   test_ = "clear() function";

   // Clearing a default constructed vector
   {
      blaze::StaticVector<int,4UL,blaze::rowVector> vec;

      clear( vec );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 0UL );

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Clearing an initialized vector
   {
      // Initialization check
      blaze::StaticVector<int,4UL,blaze::rowVector> vec{ 1, 2, 3, 4 };

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 4UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 3 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 3 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing a single element
      clear( vec[2] );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 3UL );

      if( vec[0] != 1 || vec[1] != 2 || vec[2] != 0 || vec[3] != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 1 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }

      // Clearing the vector
      clear( vec );

      checkSize    ( vec, 4UL );
      checkCapacity( vec, 4UL );
      checkNonZeros( vec, 0UL );

      if( vec[0] != 0 || vec[1] != 0 || vec[2] != 0 || vec[3] != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Clear operation failed\n"
             << " Details:\n"
             << "   Result:\n" << vec << "\n"
             << "   Expected result:\n( 0 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the StaticVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the StaticVector class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   test_ = "StaticVector swap";

   blaze::StaticVector<int,4UL,blaze::rowVector> vec1{ 1, 2, 3, 4 };
   blaze::StaticVector<int,4UL,blaze::rowVector> vec2{ 4, 3, 2, 1 };

   swap( vec1, vec2 );

   checkSize    ( vec1, 4UL );
   checkCapacity( vec1, 4UL );
   checkNonZeros( vec1, 4UL );

   if( vec1[0] != 4 || vec1[1] != 3 || vec1[2] != 2 || vec1[3] != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the first vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 4 3 2 1 )\n";
      throw std::runtime_error( oss.str() );
   }

   checkSize    ( vec2, 4UL );
   checkCapacity( vec2, 4UL );
   checkNonZeros( vec2, 4UL );

   if( vec2[0] != 1 || vec2[1] != 2 || vec2[2] != 3 || vec2[3] != 4 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the second vector failed\n"
          << " Details:\n"
          << "   Result:\n" << vec1 << "\n"
          << "   Expected result:\n( 1 2 3 4 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c isDefault() function with the StaticVector class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c isDefault() function with the StaticVector class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIsDefault()
{
   using blaze::isDefault;

   test_ = "isDefault() function";

   // isDefault with default vector
   {
      blaze::StaticVector<int,3UL,blaze::rowVector> vec;

      if( isDefault( vec ) != true ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // isDefault with non-default vector
   {
      blaze::StaticVector<int,3UL,blaze::rowVector> vec{ 0, 1, 0 };

      if( isDefault( vec[1] ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector element: " << vec[1] << "\n";
         throw std::runtime_error( oss.str() );
      }

      if( isDefault( vec ) != false ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid isDefault evaluation\n"
             << " Details:\n"
             << "   Vector:\n" << vec << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace staticvector

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
   std::cout << "   Running StaticVector class test..." << std::endl;

   try
   {
      RUN_STATICVECTOR_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during StaticVector class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
