//=================================================================================================
/*!
//  \file src/mathtest/matrices/dynamicmatrix/ClassTest1.cpp
//  \brief Source file for the DynamicMatrix class test (part 1)
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
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/CustomMatrix.h>
#include <blaze/math/DiagonalMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blaze/math/UpperMatrix.h>
#include <blaze/util/Complex.h>
#include <blaze/util/Memory.h>
#include <blaze/util/policies/Deallocate.h>
#include <blaze/util/Random.h>
#include <blazetest/mathtest/matrices/dynamicmatrix/ClassTest.h>
#include <blazetest/mathtest/RandomMaximum.h>
#include <blazetest/mathtest/RandomMinimum.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace dynamicmatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the DynamicMatrix class test.
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
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the DynamicMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the DynamicMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   //=====================================================================================
   // Row-major default constructor
   //=====================================================================================

   // Default constructor
   {
      test_ = "Row-major DynamicMatrix default constructor";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Row-major size constructor
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix size constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix size constructor (0x4)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 4UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix size constructor (3x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix size constructor (3x4)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
   }


   //=====================================================================================
   // Row-major homogeneous initialization
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix homogeneous initialization constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 0UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix homogeneous initialization constructor (0x4)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 0UL, 4UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix homogeneous initialization constructor (3x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 0UL, 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix homogeneous initialization constructor (3x4)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL, 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat, 12UL );
      checkNonZeros( mat,  0UL, 4UL );
      checkNonZeros( mat,  1UL, 4UL );
      checkNonZeros( mat,  2UL, 4UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major list initialization
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix initializer list constructor (complete list)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1, 2, 3 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major DynamicMatrix initializer list constructor (incomplete list)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat{ { 1 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array initialization
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix dynamic array initialization constructor";

      std::unique_ptr<int[]> array( new int[6] );
      array[0] = 1;
      array[1] = 2;
      array[2] = 3;
      array[3] = 4;
      array[4] = 5;
      array[5] = 6;
      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 2UL, 3UL, array.get() );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major DynamicMatrix static array initialization constructor";

      const int array[2][3] = { { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::rowMajor> mat( array );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major DynamicMatrix std::array initialization constructor";

      const std::array<std::array<int,3UL>,2UL> array{ { { 1, 2, 3 }, { 4, 5, 6 } } };
      blaze::DynamicMatrix<int,blaze::rowMajor> mat( array );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy constructor
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix copy constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 0UL, 0UL );
      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix copy constructor (0x3)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 0UL, 3UL );
      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix copy constructor (2x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 2UL, 0UL );
      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix copy constructor (2x3)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 1, 2, 3 },
                                                      { 4, 5, 6 } };

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major move constructor
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix move constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 0UL, 0UL );
      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix move constructor (0x3)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 0UL, 3UL );
      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix move constructor (2x0)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1( 2UL, 0UL );
      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Row-major DynamicMatrix copy constructor (2x3)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 1, 2, 3 },
                                                      { 4, 5, 6 } };

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix constructor
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix constructor (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix constructor (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix constructor (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix constructor (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix constructor
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      const blaze::DynamicMatrix<int,blaze::rowMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major default constructor
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix default constructor";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat;

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }


   //=====================================================================================
   // Column-major size constructor
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix size constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 0UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix size constructor (0x4)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 4UL );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix size constructor (3x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix size constructor (3x4)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 4UL );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
   }


   //=====================================================================================
   // Column-major homogeneous initialization
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix homogeneous initialization constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 0UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix homogeneous initialization constructor (0x4)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 0UL, 4UL, 2 );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix homogeneous initialization constructor (3x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 0UL, 2 );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix homogeneous initialization constructor (3x4)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 4UL, 2 );

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat, 12UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );
      checkNonZeros( mat,  2UL, 3UL );
      checkNonZeros( mat,  3UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list initialization
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix initializer list constructor (complete list)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1, 2, 3 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major DynamicMatrix initializer list constructor (incomplete list)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat{ { 1 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major array initialization
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix dynamic array initialization constructor";

      std::unique_ptr<int[]> array( new int[6] );
      array[0] = 1;
      array[1] = 2;
      array[2] = 3;
      array[3] = 4;
      array[4] = 5;
      array[5] = 6;
      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 2UL, 3UL, array.get() );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 3 || mat(0,2) != 5 ||
          mat(1,0) != 2 || mat(1,1) != 4 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 3 5 )\n( 2 4 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major DynamicMatrix static array initialization constructor";

      const int array[2][3] = { { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::columnMajor> mat( array );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major DynamicMatrix std::array initialization constructor";

      const std::array<std::array<int,3UL>,2UL> array{ { { 1, 2, 3 }, { 4, 5, 6 } } };
      blaze::DynamicMatrix<int,blaze::columnMajor> mat( array );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy constructor
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix copy constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 0UL, 0UL );
      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix copy constructor (0x3)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 0UL, 3UL );
      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix copy constructor (2x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 2UL, 0UL );
      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix copy constructor (2x3)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 1, 2, 3 },
                                                         { 4, 5, 6 } };

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major move constructor
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix move constructor (0x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 0UL, 0UL );
      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix move constructor (0x3)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 0UL, 3UL );
      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 0UL );
      checkColumns ( mat2, 3UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix move constructor (2x0)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1( 2UL, 0UL );
      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 0UL );
      checkNonZeros( mat2, 0UL );
   }

   {
      test_ = "Column-major DynamicMatrix move constructor (2x3)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 1, 2, 3 },
                                                         { 4, 5, 6 } };

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( std::move( mat1 ) );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix constructor
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix constructor (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix constructor (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix constructor (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix constructor (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix constructor
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix constructor";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      const blaze::DynamicMatrix<int,blaze::columnMajor> mat2( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Construction failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DynamicMatrix assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all assignment operators of the DynamicMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAssignment()
{
   //=====================================================================================
   // Row-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix homogeneous assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat( 3UL, 4UL );
      mat = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat, 12UL );
      checkNonZeros( mat,  0UL, 4UL );
      checkNonZeros( mat,  1UL, 4UL );
      checkNonZeros( mat,  2UL, 4UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major list assignment
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix initializer list assignment (complete list)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat;
      mat = { { 1, 2, 3 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major DynamicMatrix initializer list assignment (incomplete list)";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat;
      mat = { { 1 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 1UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major array assignment
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix static array assignment";

      const int array[2][3] = { { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::rowMajor> mat;
      mat = array;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major DynamicMatrix std::array assignment";

      const std::array<std::array<int,3UL>,2UL> array{ { { 1, 2, 3 }, { 4, 5, 6 } } };
      blaze::DynamicMatrix<int,blaze::rowMajor> mat;
      mat = array;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 3UL );
      checkNonZeros( mat, 1UL, 3UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major copy assignment
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix copy assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 1, 2, 3 },
                                                      { 4, 5, 6 } };

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major DynamicMatrix copy assignment stress test";

      using RandomMatrixType = blaze::DynamicMatrix<int,blaze::rowMajor>;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Row-major move assignment
   //=====================================================================================

   {
      test_ = "Row-major DynamicMatrix move assignment";

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1{ { 1, 2, 3 },
                                                      { 4, 5, 6 } };

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 11 }, { 12 }, { 13 }, { 14 } };

      mat2 = std::move( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::rowMajor> mat1{ { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix assignment stress test";

      using RandomMatrixType = blaze::DynamicMatrix<short,blaze::rowMajor>;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1;
      const short min( randmin );
      const short max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::columnMajor> mat1{ { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 3UL );
      checkNonZeros( mat2, 1UL, 3UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix assignment stress test";

      using RandomMatrixType = blaze::DynamicMatrix<short,blaze::columnMajor>;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1;
      const short min( randmin );
      const short max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix assignment stress test";

      using RandomMatrixType = blaze::CompressedMatrix<int,blaze::rowMajor>;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix assignment stress test";

      using RandomMatrixType = blaze::CompressedMatrix<int,blaze::columnMajor>;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major homogeneous assignment
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix homogeneous assigment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat( 3UL, 4UL );
      mat = 2;

      checkRows    ( mat,  3UL );
      checkColumns ( mat,  4UL );
      checkCapacity( mat, 12UL );
      checkNonZeros( mat, 12UL );
      checkNonZeros( mat,  0UL, 3UL );
      checkNonZeros( mat,  1UL, 3UL );
      checkNonZeros( mat,  2UL, 3UL );
      checkNonZeros( mat,  3UL, 3UL );

      if( mat(0,0) != 2 || mat(0,1) != 2 || mat(0,2) != 2 || mat(0,3) != 2 ||
          mat(1,0) != 2 || mat(1,1) != 2 || mat(1,2) != 2 || mat(1,3) != 2 ||
          mat(2,0) != 2 || mat(2,1) != 2 || mat(2,2) != 2 || mat(2,3) != 2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 2 2 2 2 )\n( 2 2 2 2 )\n( 2 2 2 2 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major list assignment
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix initializer list assignment (complete list)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat;
      mat = { { 1, 2, 3 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major DynamicMatrix initializer list assignment (incomplete list)";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat;
      mat = { { 1 }, { 4, 5, 6 } };

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 4UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 1UL );
      checkNonZeros( mat, 2UL, 1UL );

      if( mat(0,0) != 1 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 0 0 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major array assignment
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix static array initialization constructor";

      const int array[2][3] = { { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::columnMajor> mat;
      mat = array;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major DynamicMatrix std::array initialization constructor";

      const std::array<std::array<int,3UL>,2UL> array{ { { 1, 2, 3 }, { 4, 5, 6 } } };
      blaze::DynamicMatrix<int,blaze::columnMajor> mat;
      mat = array;

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 6UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 2UL );
      checkNonZeros( mat, 2UL, 2UL );

      if( mat(0,0) != 1 || mat(0,1) != 2 || mat(0,2) != 3 ||
          mat(1,0) != 4 || mat(1,1) != 5 || mat(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major copy assignment
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix copy assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 1, 2, 3 },
                                                         { 4, 5, 6 } };

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major DynamicMatrix copy assignment stress test";

      using RandomMatrixType = blaze::DynamicMatrix<int,blaze::columnMajor>;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }


   //=====================================================================================
   // Column-major move assignment
   //=====================================================================================

   {
      test_ = "Column-major DynamicMatrix move assignment";

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1{ { 1, 2, 3 },
                                                         { 4, 5, 6 } };

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 11 }, { 12 }, { 13 }, { 14 } };

      mat2 = std::move( mat1 );

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::rowMajor> mat1{ { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix assignment stress test";

      using RandomMatrixType = blaze::DynamicMatrix<short,blaze::rowMajor>;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1;
      const short min( randmin );
      const short max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::columnMajor> mat1{ { 1, 2, 3 }, { 4, 5, 6 } };
      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(0,2) = 3;
      mat1(1,0) = 4;
      mat1(1,1) = 5;
      mat1(1,2) = 6;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 6UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 3 ||
          mat2(1,0) != 4 || mat2(1,1) != 5 || mat2(1,2) != 6 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 3 )\n( 4 5 6 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix assignment stress test";

      using RandomMatrixType = blaze::DynamicMatrix<short,blaze::columnMajor>;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1;
      const short min( randmin );
      const short max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix assignment stress test";

      using RandomMatrixType = blaze::CompressedMatrix<int,blaze::rowMajor>;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL );
      mat1(0,0) = 1;
      mat1(0,1) = 2;
      mat1(1,0) = 3;
      mat1(1,2) = 4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2;
      mat2 = mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 1UL );
      checkNonZeros( mat2, 2UL, 1UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 2 || mat2(0,2) != 0 ||
          mat2(1,0) != 3 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 2 0 )\n( 3 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix assignment stress test";

      using RandomMatrixType = blaze::CompressedMatrix<int,blaze::columnMajor>;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat1;
      const int min( randmin );
      const int max( randmax );

      for( size_t i=0UL; i<100UL; ++i )
      {
         const size_t rows   ( blaze::rand<size_t>( 0UL, 10UL ) );
         const size_t columns( blaze::rand<size_t>( 0UL, 10UL ) );
         const RandomMatrixType mat2( blaze::rand<RandomMatrixType>( rows, columns, min, max ) );

         mat1 = mat2;

         if( mat1 != mat2 ) {
            std::ostringstream oss;
            oss << " Test: " << test_ << "\n"
                << " Error: Assignment failed\n"
                << " Details:\n"
                << "   Result:\n" << mat1 << "\n"
                << "   Expected result:\n" << mat2 << "\n";
            throw std::runtime_error( oss.str() );
         }
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL );
      randomize( mat2 );

      mat2 = mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DynamicMatrix addition assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the addition assignment operators of the DynamicMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testAddAssign()
{
   //=====================================================================================
   // Row-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix addition assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::rowMajor> mat1{ {  1, 2, 0 },
                                                        { -3, 0, 4 } };

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix addition assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::columnMajor> mat1{ {  1, 2, 0 },
                                                           { -3, 0, 4 } };

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix addition assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::rowMajor> mat1{ {  1, 2, 0 },
                                                        { -3, 0, 4 } };

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix addition assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::columnMajor> mat1{ {  1, 2, 0 },
                                                        { -3, 0, 4 } };

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix addition assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix addition assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix addition assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix addition assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) =  1;
      mat1(0,1) =  2;
      mat1(1,0) = -3;
      mat1(1,2) =  4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 += mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix addition assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix addition assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix addition assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 += mat1;

      if( mat1 != mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Addition assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the DynamicMatrix subtraction assignment operators.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the subtraction assignment operators of the DynamicMatrix
// class template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSubAssign()
{
   //=====================================================================================
   // Row-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix subtraction assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::rowMajor> mat1{ { -1, -2,  0 },
                                                        {  3,  0, -4 } };

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix subtraction assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::columnMajor> mat1{ { -1, -2,  0 },
                                                           {  3,  0, -4 } };

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 2UL, 3UL, 0 );
      mat2(0,1) = -2;
      mat2(0,2) =  6;
      mat2(1,0) =  5;

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Row-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2{ { 0, -2, 6 },
                                                      { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/row-major DynamicMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Row-major/column-major DynamicMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::rowMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major dense matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix subtraction assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::rowMajor> mat1{ { -1, -2,  0 },
                                                        {  3,  0, -4 } };

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::rowMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,rowMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 32UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::rowMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,rowMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix subtraction assignment (mixed type)";

      blaze::DynamicMatrix<short,blaze::columnMajor> mat1{ { -1, -2,  0 },
                                                           {  3,  0, -4 } };

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix subtraction assignment (aligned/padded)";

      using blaze::aligned;
      using blaze::padded;
      using blaze::columnMajor;

      using AlignedPadded = blaze::CustomMatrix<int,aligned,padded,columnMajor>;
      std::unique_ptr<int[],blaze::Deallocate> memory( blaze::allocate<int>( 48UL ) );
      AlignedPadded mat1( memory.get(), 2UL, 3UL, 16UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix subtraction assignment (unaligned/unpadded)";

      using blaze::unaligned;
      using blaze::unpadded;
      using blaze::columnMajor;

      using UnalignedUnpadded = blaze::CustomMatrix<int,unaligned,unpadded,columnMajor>;
      std::unique_ptr<int[]> memory( new int[7UL] );
      UnalignedUnpadded mat1( memory.get()+1UL, 2UL, 3UL );
      mat1 = 0;
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix dense matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::DynamicMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }


   //=====================================================================================
   // Column-major sparse matrix subtraction assignment
   //=====================================================================================

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::rowMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix subtraction assignment";

      blaze::CompressedMatrix<int,blaze::columnMajor> mat1( 2UL, 3UL, 4UL );
      mat1(0,0) = -1;
      mat1(0,1) = -2;
      mat1(1,0) =  3;
      mat1(1,2) = -4;

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2{ { 0, -2, 6 },
                                                         { 5,  0, 0 } };

      mat2 -= mat1;

      checkRows    ( mat2, 2UL );
      checkColumns ( mat2, 3UL );
      checkCapacity( mat2, 6UL );
      checkNonZeros( mat2, 4UL );
      checkNonZeros( mat2, 0UL, 2UL );
      checkNonZeros( mat2, 1UL, 0UL );
      checkNonZeros( mat2, 2UL, 2UL );

      if( mat2(0,0) != 1 || mat2(0,1) != 0 || mat2(0,2) != 6 ||
          mat2(1,0) != 2 || mat2(1,1) != 0 || mat2(1,2) != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat2 << "\n"
             << "   Expected result:\n( 1 0 6 )\n( 2 0 4 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix subtraction assignment (lower)";

      blaze::LowerMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix subtraction assignment (upper)";

      blaze::UpperMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/row-major DynamicMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::rowMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      test_ = "Column-major/column-major DynamicMatrix sparse matrix subtraction assignment (diagonal)";

      blaze::DiagonalMatrix< blaze::CompressedMatrix<int,blaze::columnMajor> > mat1( 3UL );
      randomize( mat1 );

      blaze::DynamicMatrix<int,blaze::columnMajor> mat2( 3UL, 3UL, 0 );

      mat2 -= mat1;

      if( mat1 != -mat2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Subtraction assignment failed\n"
             << " Details:\n"
             << "   Result:\n" << mat1 << "\n"
             << "   Expected result:\n" << mat2 << "\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************

} // namespace dynamicmatrix

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
   std::cout << "   Running DynamicMatrix class test (part 1)..." << std::endl;

   try
   {
      RUN_DYNAMICMATRIX_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during DynamicMatrix class test (part 1):\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
